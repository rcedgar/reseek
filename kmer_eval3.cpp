#include "myutils.h"
#include "scop40bench.h"
#include "cigar.h"

static void GetColPosVecs(const string &Path, char DI, uint Lo,
  vector<uint> &ColToPos, vector<uint> &PosToCol)
	{
	ColToPos.clear();
	PosToCol.clear();
	const uint ColCount = SIZE(Path);
	for (uint i = 0; i < Lo; ++i)
		PosToCol.push_back(UINT_MAX);
	uint Pos = Lo;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == DI)
			{
			ColToPos.push_back(Pos);
			PosToCol.push_back(Col);
			}
		else
			ColToPos.push_back(UINT_MAX);
		++Pos;
		}
	}

static bool HasSeedOnPath(const vector<uint> &SeedPosQs,
  const vector<uint> &SeedPosRs, uint LoQ, uint LoR, const string &Path)
	{
	vector<uint> ColToPosQ;
	vector<uint> ColToPosR;
	vector<uint> PosToColQ;
	vector<uint> PosToColR;
	GetColPosVecs(Path, 'D', LoQ, ColToPosQ, PosToColQ);
	GetColPosVecs(Path, 'I', LoR, ColToPosR, PosToColR);
	for (uint i = 0; i < SIZE(SeedPosQs); ++i)
		{
		uint PosQ = SeedPosQs[i];
		uint PosR = SeedPosRs[i];
	// Local alignment, ok to overflow here
		if (PosQ >= SIZE(PosToColQ))
			continue;
		if (PosR >= SIZE(PosToColR))
			continue;
		uint ColQ = PosToColQ[PosQ];
		uint ColR = PosToColR[PosR];
		if (ColQ == ColR && ColQ != UINT_MAX)
			return true;
		}
	return false;
	}

void cmd_kmer_eval3()
	{
	const string &TsvFN = g_Arg1;
	const string &CalFN = opt_input;

	DSSParams Params;
	Params.SetFromCmdLine();
	Params.m_MAXNQNR = 32;
	Params.m_MAXFX = 8;
	Params.m_MinDiagSpacing = 8;
	Params.m_MaxDiagSpacing = 32;
	Params.m_MinKmerScore = -99; // 0.01f;
	Params.m_PatternStr = "101000011";
	Params.m_X = 4;

	DSSAligner DA;
	DA.m_Params = &Params;

	SCOP40Bench SB;
	SB.Setup(Params);
	SB.ReadChains(CalFN);

	vector<double> Evalues;
	vector<bool> Ts;
	vector<uint> Dom1s;
	vector<uint> Dom2s;
	vector<uint> Lo1s;
	vector<uint> Lo2s;
	vector<string> Paths;

	FILE *f = OpenStdioFile(TsvFN);
	string Line;
	vector<string> Fields;
	uint NT = 0;
	uint NF = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 5);
		asserta(SIZE(Fields[0]) == 1);
		char TF = Fields[0][0];
		double Evalue = StrToFloat(Fields[1]);
		const string &Dom1 = Fields[2];
		const string &Dom2 = Fields[3];
		const string &CIGAR = Fields[4];
	
		string Path;
		uint Lo1, Lo2;
		LocalCIGARToPath(CIGAR, Path, Lo1, Lo2, false);

		uint DomIdx1 = SB.GetDomIdx(Dom1);
		uint DomIdx2 = SB.GetDomIdx(Dom2);

		Ts.push_back(TF == 'T');
		if (TF == 'T')
			++NT;
		else if (TF == 'F')
			++NF;
		else
			asserta(false);
		Evalues.push_back(Evalue);
		Dom1s.push_back(DomIdx1);
		Dom2s.push_back(DomIdx2);
		Lo1s.push_back(Lo1);
		Lo2s.push_back(Lo2);
		Paths.push_back(Path);

		uint PathL1, PathL2;
		PathToLs(Path, PathL1, PathL2);

		const PDBChain &Chain1 = SB.GetChainByDomIdx(DomIdx1);
		const PDBChain &Chain2 = SB.GetChainByDomIdx(DomIdx2);
		uint L1 = Chain1.GetSeqLength();
		uint L2 = Chain2.GetSeqLength();
		asserta(PathL1 <= L1);
		asserta(PathL2 <= L2);

		const vector<vector<byte> > &Profile1 = SB.GetProfileByDomIdx(DomIdx1);
		const vector<vector<byte> > &Profile2 = SB.GetProfileByDomIdx(DomIdx2);

		float Evalue2 = DA.GetEvaluePath(Chain1, Chain2,
		  Profile1, Profile2, Lo1, Lo2, Path);

		asserta(feq(Evalue, Evalue2, Evalue/100));
		}
	CloseStdioFile(f);

	const uint N = SIZE(Dom1s);
	uint PairsWithDiagT = 0;
	uint PairsWithDiagF = 0;
	uint Counter = 0;
	uint TPairsWithPveU = 0;
	uint TPairsWithSOP = 0;
	uint FPairsWithPveU = 0;
	const uint UBINS = 32;
	vector<uint> UToCountT(UBINS+1);
	vector<uint> UToCountF(UBINS+1);
	for (uint i = 0; i < N; ++i)
		{
		uint DomIdxQ = Dom1s[i];
		uint DomIdxR = Dom2s[i];
		uint Lo1 = Lo1s[i];
		uint Lo2 = Lo2s[i];
		bool T = Ts[i];
		const string &Path = Paths[i];

		const PDBChain &ChainQ = SB.GetChainByDomIdx(DomIdxQ);
		const PDBChain &ChainR = SB.GetChainByDomIdx(DomIdxR);

		const vector<vector<byte> > &ProfileQ = SB.m_Profiles[DomIdxQ];
		const vector<vector<byte> > &ProfileR = SB.m_Profiles[DomIdxR];

		const vector<uint> &KmersQ = SB.m_KmersVec[DomIdxQ];
		const vector<uint> &KmersR = SB.m_KmersVec[DomIdxR];

		vector<uint> SeedPosQs;
		vector<uint> SeedPosRs;
		DA.GetSeeds(ProfileQ, ProfileR, KmersQ, KmersR, SeedPosQs, SeedPosRs);
		if (T)
			{
			bool HasSOP = HasSeedOnPath(SeedPosQs, SeedPosRs, Lo1, Lo2, Path);
			if (HasSOP)
				++TPairsWithSOP;
			}

		uint U = SIZE(SeedPosQs); // DA.GetU(KmersQ, KmersR);
		if (U > 0)
			{
			if (T)
				++TPairsWithPveU;
			else
				++FPairsWithPveU;
			}
		if (U > UBINS)
			U = UBINS;
		if (T)
			UToCountT[U] += 1;
		else
			UToCountF[U] += 1;
		}
	ProgressLog("Pattern %s T %u(%.1f%%), F %u(%.1f%%) TSOP %u(%.1f%%)\n",
	  Params.m_PatternStr.c_str(),
	  TPairsWithPveU, GetPct(TPairsWithPveU, NT),
	  FPairsWithPveU, GetPct(FPairsWithPveU, NF),
	  TPairsWithSOP, GetPct(TPairsWithSOP, NT));

	for (uint U = 0; U <= UBINS; ++U)
		{
		Log("%2u  %5u  %5u\n", U, UToCountT[U], UToCountF[U]);
		}
	}

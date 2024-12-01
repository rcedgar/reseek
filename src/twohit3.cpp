#include "myutils.h"
#include "seqdb.h"
#include "twohitdiag.h"
#include "mudex.h"
#include "mermx.h"
#include "diaghsp.h"
#include "diag.h"
#include "alpha.h"
#include "scop40bench.h"
#include "dssaligner.h"

#define MAP_CHECK	0	// Verify by keeping count map

#if MAP_CHECK
struct DiagInfo
	{
	uint SeqIdxQ = UINT_MAX;
	uint SeqIdxT = UINT_MAX;
	uint16_t Diag = UINT16_MAX;
	bool operator<(const DiagInfo &rhs) const
		{
		if (SeqIdxQ < rhs.SeqIdxQ)
			return true;
		if (SeqIdxQ > rhs.SeqIdxQ)
			return false;
		if (SeqIdxT < rhs.SeqIdxT)
			return true;
		if (SeqIdxT > rhs.SeqIdxT)
			return false;
		return Diag < rhs.Diag;
		}
	};
#endif

void StrToMuLetters(const string &StrSeq, byte *Letters);

void cmd_twohit3()
	{
	asserta(optset_hsminscore);
	asserta(optset_mindiagscore);
	asserta(optset_ref);

	vector<PDBChain *> Chains;
	ReadChains(opt_ref, Chains);
	const uint ChainCount = SIZE(Chains);
	uint DBSize = ChainCount;
	if (optset_dbsize)
		DBSize = uint(opt_dbsize);
	ProgressLog("dbsize=%u\n", DBSize);

	DSSParams Params;
	Params.SetFromCmdLine(DBSize);

	DSS D;
	D.SetParams(Params);

	DSSAligner DA;
	DA.m_Params = &Params;

	vector<vector<vector<byte> > *> Profiles(ChainCount);
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "Profiles");
		D.Init(*Chains[i]);
		vector<vector<byte> > *Profile = new vector<vector<byte> >;
		D.GetProfile(*Profile);
		Profiles[i] = Profile;
		}

	const int HSMinScore = opt_hsminscore;
	const int MinDiagScore = opt_mindiagscore;
	const string &RefFN = opt_ref; // "c:/int/reseek_bench/alns/reseek_fast.tsv";

	bool Self = false;
	FILE *fOut = CreateStdioFile(opt_output);

#if MAP_CHECK
	map<DiagInfo, uint> DiagInfoToCount;
#endif

	SeqDB Input;
	Progress("Read Mu fasta...");
	Input.FromFasta(g_Arg1);
	Input.SetLabelToIndex();
	Progress(" done.\n");
	const uint SeqCount = Input.GetSeqCount();

	byte **MuSeqs = myalloc(byte *, SeqCount);
	vector<vector<byte> > MuLettersVec(SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "Mu letters");
		const string &Label = Input.GetLabel(SeqIdx);
		const string &ChainLabel = Chains[SeqIdx]->m_Label;
		asserta(Label == ChainLabel);
		uint L = Input.GetSeqLength(SeqIdx);
		byte *MuSeq = myalloc(byte, L);
		StrToMuLetters(Input.GetSeq(SeqIdx), MuSeq);
		MuSeqs[SeqIdx] = MuSeq;

		vector<byte> &MuLetters = MuLettersVec[SeqIdx];
		MuLetters.reserve(L);
		for (uint i = 0; i < L; ++i)
			MuLetters.push_back(MuSeq[i]);
		}

	MuDex MD;
	MD.FromSeqDB(Input);
	const uint k = MD.m_k;
	const uint DictSize = MD.m_DictSize;

	MerMx MM;
	extern const short * const *Mu_S_ij_short;
	MM.Init(Mu_S_ij_short, 5, 36, 2);
	asserta(MM.m_AS_pow[5] == DictSize);

	DiagHSP DH;

// Buffer for high-scoring k-mers
	uint *HSKmers = myalloc(uint, DictSize);

	int *SeqIdxTToBestDiagScore = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiag = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiagLo = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiagLen = myalloc(int, SeqCount);
	uint *SeqIdxTs = myalloc(uint, SeqCount);
	uint HitCount = 0;

	zero_array(SeqIdxTToBestDiagScore, SeqCount);
	zero_array(SeqIdxTToBestDiag, SeqCount);
	zero_array(SeqIdxTToBestDiagLo, SeqCount);
	zero_array(SeqIdxTToBestDiagLen, SeqCount);
	zero_array(SeqIdxTs, SeqCount);

	uint NH2 = 0;
	uint NT2 = 0;
	uint NT3 = 0;
	TwoHitDiag TH;
	for (uint SeqIdxQ = 0; SeqIdxQ < SeqCount; ++SeqIdxQ)
		{
		double FractDone = double(SeqIdxQ)/SeqCount;
		double IdealNT = FractDone*173058;
		double PctIdeal = NT3*100.0/(IdealNT + 1);
		ProgressStep(SeqIdxQ, SeqCount, "Searching %u %u (%.1f, %.1f%%)",
					 NT2, NT3, IdealNT, PctIdeal);
		const PDBChain &ChainQ = *Chains[SeqIdxQ];
		const vector<vector<byte> > *ptrProfileQ = Profiles[SeqIdxQ];
		vector<byte> *ptrMuLettersQ = &MuLettersVec[SeqIdxQ];
		DA.SetQuery(ChainQ, ptrProfileQ, 0, ptrMuLettersQ, 0);

		const char *SeqQ = Input.GetSeq(SeqIdxQ).c_str();
		const string &LabelQ = Input.GetLabel(SeqIdxQ);
		uint LQ32 = Input.GetSeqLength(SeqIdxQ);
		asserta(LQ32 < UINT16_MAX);
		uint16_t LQ = uint16_t(LQ32);

	// Initialize k-mer scan of Query
		uint Kmer = 0;
		for (uint SeqPosQ = 0; SeqPosQ < k-1; ++SeqPosQ)
			{
			byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
			assert(Letter < 36);
			Kmer = Kmer*36 + Letter;
			}

	// k-mer scan of Query
		for (uint SeqPosQ = k-1; SeqPosQ < LQ; ++SeqPosQ)
			{
			byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
			assert(Letter < 36);
			Kmer = Kmer*36 + Letter;
			Kmer %= DictSize;
		// Construct high-scoring neighborhood of current k-mer (Kmer)
			const uint HSKmerCount =
				MM.GetHighScoring5mers(Kmer, HSMinScore, HSKmers);

			for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
				{
				uint HSKmer = HSKmers[HSKmerIdx];
			// Look up HSKmer in MuDex (Mu k-mer index of db)
				uint RowSize = MD.GetRowSize(HSKmer);
				uint DataOffset = MD.GetRowStart(HSKmer);
				for (uint ColIdx = 0; ColIdx < RowSize; ++ColIdx)
					{
					uint32_t SeqIdxT;
					uint16_t SeqPosT;
					MD.Get(DataOffset++, SeqIdxT, SeqPosT);
					if (SeqIdxT == SeqIdxQ)
						continue;
					if (Self && SeqIdxT < SeqIdxQ)
						continue;
					uint LT32 = Input.GetSeqLength(SeqIdxT);
					asserta(LT32 < UINT16_MAX);
					uint16_t LT = uint16_t(LT32);
					diag dg(LQ, LT);
					uint16_t Diag = dg.getd(SeqPosQ, SeqPosT);
					TH.Add(SeqIdxT, Diag);
					}
				}
			}
		TH.ClearDupes();
		TH.SetDupes();
#if DEBUG
		TH.Validate(SeqCount, INT16_MAX);
#endif
		vector<pair<uint32_t, uint16_t> > SeqIdxDiagPairs;
		uint DupeCount = TH.m_DupeCount;
		const byte *MuSeqQ = MuSeqs[SeqIdxQ];
		DH.SetQ(MuSeqQ, LQ);
		assert(HitCount == 0);
		for (uint i = 0; i < DupeCount; ++i)
			{
			uint32_t SeqIdxT = TH.m_DupeSeqIdxs[i];
			uint16_t Diag = TH.m_DupeDiags[i];
			const byte *MuSeqT = MuSeqs[SeqIdxT];
			const uint LT = Input.GetSeqLength(SeqIdxT);
			DH.SetT(MuSeqT, LT);
			int Lo, Len;
			int DiagScore = DH.Search(Diag, Lo, Len);
			if (DiagScore >= MinDiagScore)
				{
				int BestDiagScoreT = SeqIdxTToBestDiagScore[SeqIdxT];
				if (BestDiagScoreT == 0)
					{
					SeqIdxTs[HitCount++] = SeqIdxT;
					SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
					SeqIdxTToBestDiag[SeqIdxT] = Diag;
					SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
					SeqIdxTToBestDiagLen[SeqIdxT] = Len;
					}
				else if (DiagScore > SeqIdxTToBestDiagScore[SeqIdxT])
					{
					SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
					SeqIdxTToBestDiag[SeqIdxT] = Diag;
					SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
					SeqIdxTToBestDiagLen[SeqIdxT] = Len;
					}
				}
			}

		if (HitCount > 0)
			{
			for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
				{
				uint SeqIdxT = SeqIdxTs[HitIdx];
				const string &LabelT = Input.GetLabel(SeqIdxT);
				bool IsT = SCOP40Bench::IsTP_SF(LabelQ, LabelT);
				++NH2;
				if (IsT)
					++NT2;
				//int BestDiagScore = SeqIdxTToBestDiagScore[SeqIdxT];
				//int BestDiag = SeqIdxTToBestDiag[SeqIdxT];
				//int BestDiagLo = SeqIdxTToBestDiagLo[SeqIdxT];
				//int BestDiagLen = SeqIdxTToBestDiagLen[SeqIdxT];

				const PDBChain *ChainT = Chains[SeqIdxT];
				const vector<vector<byte> > *ptrProfileT = Profiles[SeqIdxT];
				vector<byte> *ptrMuLettersT = &MuLettersVec[SeqIdxT];
				DA.SetTarget(*ChainT, ptrProfileT, 0, ptrMuLettersT, 0);
				DA.AlignQueryTarget();
				float Evalue = DA.m_EvalueA;
				if (IsT && Evalue <= 10)
					++NT3;
				if (fOut != 0)
					fprintf(fOut, "%s\t%s\t%.3g\n",
							LabelQ.c_str(), LabelT.c_str(), Evalue);
				}
			}

		for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
			{
			uint SeqIdxT = SeqIdxTs[HitIdx];
			SeqIdxTToBestDiagScore[SeqIdxT] = 0;
			}
#if DEBUG
		{
		for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
			{
			assert(SeqIdxTToBestDiagScore[SeqIdx] == 0);
			}
		}
#endif
		HitCount = 0;
		TH.Reset();
		}
	CloseStdioFile(fOut);

	double PctIdeal = NT3*100.0/173058;

	ProgressLog("NT2=%u", NT2);
	ProgressLog("\tNT3=%u", NT3);
	ProgressLog("\tIdeal=%.1f%%", PctIdeal);
	ProgressLog("\tNH=%u", NH2);
	ProgressLog("\ths=%u", opt_hsminscore);
	ProgressLog("\tdg=%u", opt_mindiagscore);
	}

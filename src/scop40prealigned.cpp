#include "myutils.h"
#include "scop40bench.h"
#include "seqdb.h"
#include "alpha.h"
#include "cigar.h"

void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);

static FILE *s_fIn;
static FILE *s_fOut;
static mutex s_InputLock;
static mutex s_OutpuLock;

static void GetPathFromRows(const string &Row1, const string &Row2,
							string &Path)
	{
	Path.clear();
	const uint ColCount = SIZE(Row1);
	asserta(SIZE(Row2) == ColCount);
	uint Pos1 = 0;
	uint Pos2 = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c1 = Row1[Col];
		char c2 = Row2[Col];
		bool gap1 = isgap(c1);
		bool gap2 = isgap(c2);
		if (!gap1 && !gap2)
			Path += 'M';
		else if (!gap1 && gap2)
			Path += 'D';
		else if (gap1 && !gap2)
			Path += 'I';
		else
			asserta(false);
		}
	}

void SCOP40Bench::StaticThreadBodyPreAligned(uint ThreadIndex, SCOP40Bench *ptrSB)
	{
	ptrSB->ThreadBodyPrealigned(ThreadIndex);
	}

void SCOP40Bench::ThreadBodyPrealigned(uint ThreadIndex)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	DSSAligner DA;
	DA.SetParams(*m_Params);

	string Line;
	vector<string> Fields;
	for (;;)
		{
		bool Ok = false;
		s_InputLock.lock();
		Ok = ReadLineStdioFile(s_fIn, Line);
		s_InputLock.unlock();
		if (!Ok)
			return;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 4);

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		float OldDPScore = StrToFloatf(Fields[2]);
		const string &CIGAR = Fields[3];

		string Dom1, Dom2;
		GetDomFromLabel(Label1, Dom1);
		GetDomFromLabel(Label2, Dom2);

		bool IsTP = IsTP_SF(Label1, Label2);

		map<string, uint>::const_iterator iter1 = m_DomToChainIdx.find(Dom1);
		map<string, uint>::const_iterator iter2 = m_DomToChainIdx.find(Dom2);
		asserta(iter1 != m_DomToChainIdx.end());
		asserta(iter2 != m_DomToChainIdx.end());

		uint ChainIdx1 = iter1->second;
		uint ChainIdx2 = iter2->second;

		const PDBChain &Chain1 = *m_DBChains[ChainIdx1];
		const PDBChain &Chain2 = *m_DBChains[ChainIdx2];
		uint L1 = Chain1.GetSeqLength();
		uint L2 = Chain2.GetSeqLength();

		const vector<vector<byte> > &Profile1 = *m_DBProfiles[ChainIdx1];
		const vector<vector<byte> > &Profile2 = *m_DBProfiles[ChainIdx2];

		string Path;
		CIGARToPath(CIGAR, Path);

		uint NM, ND, NI;
		GetPathCounts(Path, NM, ND, NI);
		uint Sum1 = NM + ND;
		uint Sum2 = NM + NI;
		asserta(Sum1 == L1);
		asserta(Sum2 == L2);

		float NewDPScore = DA.GetDPScoreGivenPath(Profile1, Profile2, Path);
		s_OutpuLock.lock();
		StoreScore(ChainIdx1, ChainIdx2, NewDPScore);
		s_OutpuLock.unlock();

		if (s_fOut != 0)
			{
			s_OutpuLock.lock();
			fprintf(s_fOut, "%s\t%s\t%.3g\t%.3g\t%c\n",
					Dom1.c_str(),
					Dom2.c_str(),
					OldDPScore,
					NewDPScore,
					tof(IsTP));
			s_OutpuLock.unlock();
			}
		}
	}

void SCOP40Bench::RunPrealigned(const string &TsvFN)
	{
	s_fIn = OpenStdioFile(TsvFN);
	const uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBodyPreAligned, ThreadIndex, this);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	CloseStdioFile(s_fIn);
	}

void cmd_scop40bench_pre()
	{
	const string CalFN = g_Arg1;
	s_fOut = CreateStdioFile(opt(output));

	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption, SCOP40_DBSIZE);

	SCOP40Bench SB;
	SB.m_Params = &Params;
	SB.LoadDB(CalFN);

	asserta(SB.m_Params == &Params);
	Params.m_DBSize = (float) SB.GetDBSize();

	SB.Setup();
	
	float MaxFPR = 0.005f;
	if (optset_maxfpr)
		MaxFPR = (float) opt(maxfpr);

	SB.m_QuerySelf = true;
	SB.m_ScoresAreEvalues = true;
	if (opt(scores_are_not_evalues))
		SB.m_ScoresAreEvalues = false;
	SB.RunPrealigned(opt(input2));
	CloseStdioFile(s_fOut);
	s_fOut = 0;

	ProgressLog("%u / %u mu filter discards\n",
				DSSAligner::m_MuFilterDiscardCount.load(),
				DSSAligner::m_MuFilterInputCount.load());
	SB.WriteOutput();
	SB.WriteBit(opt(savebit));
	if (optset_sens1fp_report)
		{
		FILE *f = CreateStdioFile(opt(sens1fp_report));
		SB.WriteSens1FPReport(f);
		CloseStdioFile(f);
		}
	}

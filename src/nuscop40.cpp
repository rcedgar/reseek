#include "myutils.h"
#include "pdbchain.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "triangle.h"
#include "xdpmem.h"
#include "sort.h"

static string s_FeatureName;
static uint s_AS;
static vector<vector<float> > s_ScoreMx;

static uint s_PairIndex = UINT_MAX;
static uint s_PairCount = UINT_MAX;
static uint s_NextSeqIdx1 = UINT_MAX;
static uint s_NextSeqIdx2 = UINT_MAX;
static vector<vector<byte> > s_ByteSeqs;
static vector<string> s_Labels;

static vector<uint> s_DomIdxToSeqIdx;
static vector<uint> s_SeqIdxToDomIdx;
static vector<uint> s_DomIdx1s;
static vector<uint> s_DomIdx2s;
static vector<float> s_Scores;
static mutex s_HitsLock;

static uint GetSeqCount()
	{
	return SIZE(s_ByteSeqs);
	}

static byte GetLetter(const byte *Seq)
	{
	char tmp[3];
	tmp[0] = Seq[0];
	tmp[1] = Seq[1];
	tmp[2] = 0;
	char *endptr;
	uint Letter32 = strtoul(tmp, &endptr, 16);
	byte Letter = byte(Letter32);
	asserta(uint(Letter) == Letter32);
	return Letter;
	}

static void GetByteSeq(SeqInfo *SI, vector<byte> &ByteSeq)
	{
	ByteSeq.clear();
	const uint L2 = SI->m_L;
	const byte *Seq = SI->m_Seq;
	asserta(L2%2 == 0);
	const uint L = L2/2;
	ByteSeq.reserve(L);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		byte Letter = GetLetter(Seq + 2*Pos);
		ByteSeq.push_back(Letter);
		}
	}

static void ReadSubstMx(const string &TsvFN, string &FeatureName, uint &AS,
	vector<vector<float> > &Mx)
	{
	FeatureName.clear();
	AS = UINT_MAX;
	Mx.clear();
	FILE *f = OpenStdioFile(TsvFN);
	string Line;
	vector<string> Fields;
	bool FoundScoreMx = false;
	while (ReadLineStdioFile(f, Line))
		{
		if (Line.empty() || Line[0] == '#')
			continue;
		Split(Line, Fields, '\t');
		if (Fields[0] == "scoremx")
			{
			FoundScoreMx = true;
			break;
			}
		}
	asserta(FoundScoreMx);
	asserta(SIZE(Fields) == 3);
	FeatureName = Fields[1];
	AS = StrToUint(Fields[2]);
	asserta(AS >= 3 && AS <= 256);
	Mx.resize(AS);
	for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == AS+1);
		asserta(StrToUint(Fields[0]) == Letter1);
		vector<float> &Row = Mx[Letter1];
		for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
			Row.push_back(StrToFloatf(Fields[Letter2+1]));
		}
	CloseStdioFile(f);
	}

static bool GetNextPairSelf(uint &SeqIdx1, uint &SeqIdx2)
	{
	SeqIdx1 = UINT_MAX;
	SeqIdx2 = UINT_MAX;
	uint MyIndex = s_PairIndex++;
	if (MyIndex >= s_PairCount)
		return false;
	if (MyIndex == 0 || MyIndex+1 == s_PairCount || MyIndex%1000 == 0)
		ProgressStep(MyIndex, s_PairCount, "Aligning");
	triangle_k_to_ij(MyIndex, GetSeqCount(), SeqIdx1, SeqIdx2);
	return true;
	}

static float Align(const vector<byte> &SeqA, const vector<byte> &SeqB)
	{
	float SWFast_SubstMx(XDPMem &Mem,
		const byte *A, uint LA,
		const byte *B, uint LB,
		const vector<vector<float> > &SubstMx,
		float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
		string &Path);

	const uint LA = SIZE(SeqA);
	const uint LB = SIZE(SeqB);
	if (LA < 10 || LB < 10)
		return 0;

	XDPMem Mem;
	const byte *A = SeqA.data();
	const byte *B = SeqB.data();
	string Path;
	float Open = -3;
	float Ext = -1;
	uint Loi, Loj, Leni, Lenj;
	float Score = SWFast_SubstMx(Mem, A, LA, B, LB, s_ScoreMx,
		Open, Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

static void ThreadBody(uint ThreadIndex)
	{
	vector<uint> DomIdx1s;
	vector<uint> DomIdx2s;
	vector<float> Scores;
	for (;;)
		{
		uint SeqIdx1, SeqIdx2;
		bool Ok = GetNextPairSelf(SeqIdx1, SeqIdx2);
		if (!Ok)
			break;

		if (SeqIdx1 == SeqIdx2)
			continue;

		const vector<byte> &Seq1 = s_ByteSeqs[SeqIdx1];
		const vector<byte> &Seq2 = s_ByteSeqs[SeqIdx2];
		float Score = Align(Seq1, Seq2);
		uint DomIdx1 = s_SeqIdxToDomIdx[SeqIdx1];
		uint DomIdx2 = s_SeqIdxToDomIdx[SeqIdx2];
		DomIdx1s.push_back(SeqIdx1);
		DomIdx2s.push_back(SeqIdx2);
		Scores.push_back(Score);
		}

	s_HitsLock.lock();
	s_DomIdx1s.insert(s_DomIdx1s.end(), DomIdx1s.begin(), DomIdx1s.end());
	s_DomIdx2s.insert(s_DomIdx2s.end(), DomIdx2s.begin(), DomIdx2s.end());
	s_Scores.insert(s_Scores.end(), Scores.begin(), Scores.end());
	s_HitsLock.unlock();
	}

static void Search()
	{
	s_PairIndex = UINT_MAX;
	s_PairCount = UINT_MAX;
	s_NextSeqIdx1 = UINT_MAX;
	s_NextSeqIdx2 = UINT_MAX;

	uint SeqCount = GetSeqCount();
	s_PairIndex = 0;
	s_PairCount = SeqCount + (SeqCount*(SeqCount-1))/2;
	s_NextSeqIdx1 = 0;
	s_NextSeqIdx2 = 0;

	const uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}

void cmd_nuscop40()
	{
	const string &FastaFN = g_Arg1;
	const string &TsvFN = opt(input);

	ReadSubstMx(TsvFN, s_FeatureName, s_AS, s_ScoreMx);

	FASTASeqSource FSS;
	FSS.Open(FastaFN);
	FSS.m_AllowDigits = true;
	ObjMgr OM;
	SeqInfo *SI = OM.GetSeqInfo();
	s_ByteSeqs.resize(12000);
	s_Labels.clear();
	uint SeqCount = 0;
	while (FSS.GetNext(SI))
		{
		vector<byte> &ByteSeq = s_ByteSeqs[SeqCount++];
		GetByteSeq(SI, ByteSeq);
		s_Labels.push_back(string(SI->m_Label));
		}
	s_ByteSeqs.resize(SeqCount);
	ProgressLog("%u hex seqs\n", SeqCount);
	asserta(SIZE(s_ByteSeqs) == SIZE(s_Labels));
	Search();
	const uint HitCount = SIZE(s_Scores);
	ProgressLog("%u hits\n", HitCount);
	asserta(SIZE(s_DomIdx1s) == HitCount);
	asserta(SIZE(s_DomIdx2s) == HitCount);

	FILE *fOut = CreateStdioFile(opt(output));
	vector<uint> Order(HitCount);
	QuickSortOrderDesc(s_Scores.data(), HitCount, Order.data());
	for (uint k = 0; k < HitCount; ++k)
		{
		ProgressStep(k, HitCount, "Writing hits");
		uint i = Order[k];
		uint DomIdx1 = s_DomIdx1s[i];
		uint DomIdx2 = s_DomIdx2s[i];
		uint SeqIdx1 = s_DomIdxToSeqIdx[DomIdx1];
		uint SeqIdx2 = s_DomIdxToSeqIdx[DomIdx2];
		fprintf(fOut, "%s", s_Labels[SeqIdx1].c_str());
		fprintf(fOut, "\t%s", s_Labels[SeqIdx2].c_str());
		fprintf(fOut, "\t%.4g", s_Scores[i]);
		fprintf(fOut, "\n");
		}
	CloseStdioFile(fOut);
	}

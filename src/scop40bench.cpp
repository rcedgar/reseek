#include "myutils.h"
#include "scop40bench.h"
#include "mx.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "cigar.h"
#include "timing.h"
#include <algorithm>
#include <random>
#include <thread>
#include <set>

uint GetMatchColCount(const string &Path);
uint GetU(const vector<uint> &KmersQ, const vector<uint> &KmersR);
uint GetUBits(const vector<uint> &KmerBitsQ, const vector<uint> &KmerBitsR);

void SCOP40Bench::GetDomFromLabel(const string &Label, string &Dom)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	Dom = Fields[0];
	}

void SCOP40Bench::ParseScopLabel(const string &Label, string &Dom,
 string &Cls, string &Fold, string &SF, string &Fmy, bool MissingOk)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	if (SIZE(Fields) == 1)
		{
		if (!MissingOk)
			Die("ParseScopLabel, SCOP id missing >%s\n", Label.c_str());
		Dom = Fields[0];
		Cls = "-";
		Fold = "-";
		SF = "-";
		Fmy = "-";
		}
	else if (SIZE(Fields) == 2)
		{
		Dom = Fields[0];
		const string &ScopId = Fields[1];
	// Scop identifier looks like a.1.2.3 
	//     a    1           2      3
	// class.fold.superfamily.family
		vector<string> Fields2;
		Split(ScopId, Fields2, '.');
		if (SIZE(Fields2) != 4)
			Die("ParseScopLabel, bad SCOP id >%s\n", Label.c_str());
		Cls = Fields2[0];
		Fold = Fields2[0] + "." + Fields2[1];
		SF = Fields2[0] + "." + Fields2[1] + "." + Fields2[2];
		Fmy = Fields2[0] + "." + Fields2[1] + "." +  Fields2[2] + "." + Fields2[3];
		}
	else
		Die("ParseScopLabel, bad format >%s\n", Label.c_str());
	}

void SCOP40Bench::SetSFSizes()
	{
	GetSFSizes(m_SFSizes);
	}

uint SCOP40Bench::GetDomIdx(const string &Dom_or_DomSlashId,
  bool FailOnErr) const
	{
	vector<string> Fields;
	Split(Dom_or_DomSlashId, Fields, '/');
	const string &Dom = Fields[0];
	map<string, uint>::const_iterator iter = m_DomToIdx.find(Dom);
	if (iter == m_DomToIdx.end())
		{
		if (FailOnErr)
			Die("GetDomIdx(%s)", Dom_or_DomSlashId.c_str());
		return UINT_MAX;
		}
	uint DomIdx = iter->second;
	return DomIdx;
	}

const vector<vector<byte> > &SCOP40Bench::GetProfileByDomIdx(uint DomIdx) const
	{
	asserta(DomIdx < SIZE(m_Profiles));
	asserta(DomIdx < SIZE(m_DomIdxToChainIdx));
	uint ChainIdx = m_DomIdxToChainIdx[DomIdx];
	asserta(ChainIdx < SIZE(m_Profiles));
	return m_Profiles[ChainIdx];
	}

const PDBChain &SCOP40Bench::GetChainByDomIdx(uint DomIdx) const
	{
	asserta(DomIdx < SIZE(m_DomIdxToChainIdx));
	uint ChainIdx = m_DomIdxToChainIdx[DomIdx];
	asserta(ChainIdx < SIZE(m_Chains));
	return *m_Chains[ChainIdx];
	}


bool SCOP40Bench::KeepScore(float Score) const
	{
	if (optset_evalue && Score > opt_evalue)
		return false;
	if (m_ScoresAreEvalues)
		return Score >= 0 && Score < 100;
	else
		return Score > 0;
	}

void SCOP40Bench::OnSetup()
	{
	asserta(m_QuerySelf);
	asserta(m_ChainCount == m_QueryChainCount);
	asserta(m_DBChainCount == 0);
	asserta(optset_benchmode);
	m_Mode = string(opt_benchmode);
	asserta(m_Mode == "family" || m_Mode == "fold" || m_Mode == "ignore");
	m_ScoresAreEvalues = true;
	BuildDomSFIndexesFromQueryChainLabels();
	SetSFSizes();
	}

void SCOP40Bench::BuildDomSFIndexesFromQueryChainLabels()
	{
	m_DomIdxToChainIdx.clear();
	m_DomIdxToChainIdx.resize(m_ChainCount, UINT_MAX);
	for (uint ChainIndex = 0; ChainIndex < m_QueryChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, m_ChainCount, "Index chain labels");
		const PDBChain &Chain = *m_Chains[ChainIndex];

		string Dom;
		string Cls;
		string SF;
		string Fmy;
		string Fold;
		ParseScopLabel(Chain.m_Label, Dom, Cls, Fold, SF, Fmy);

		uint SFIdx = UINT_MAX;
		if (m_SFToIdx.find(SF) == m_SFToIdx.end())
			{
			SFIdx = SIZE(m_SFs);
			m_SFs.push_back(SF);
			m_SFToIdx[SF] = SFIdx;
			}
		else
			SFIdx = m_SFToIdx[SF];

		uint FoldIdx = UINT_MAX;
		if (m_FoldToIdx.find(Fold) == m_FoldToIdx.end())
			{
			FoldIdx = SIZE(m_Folds);
			m_Folds.push_back(Fold);
			m_FoldToIdx[Fold] = FoldIdx;
			}
		else
			FoldIdx = m_FoldToIdx[Fold];
 
		uint DomIdx = UINT_MAX;
		if (m_DomToIdx.find(Dom) != m_DomToIdx.end())
			{
			Die("Duplicate dom >%s", Chain.m_Label.c_str());
			DomIdx = m_DomToIdx[Dom];
			}
		DomIdx = SIZE(m_Doms);
		m_Doms.push_back(Dom + "/" + SF);
		m_DomToIdx[Dom] = DomIdx;
		m_DomIdxToSFIdx.push_back(SFIdx);
		m_DomIdxToFoldIdx.push_back(FoldIdx);

		m_DomIdxs.push_back(DomIdx);
		m_DomToChainIdx[Dom] = ChainIndex;
		m_DomIdxToChainIdx[DomIdx] = ChainIndex;
		}

	asserta(SIZE(m_DomIdxs) == m_QueryChainCount);
	asserta(SIZE(m_DomIdxToSFIdx) == m_QueryChainCount);
	asserta(SIZE(m_DomIdxToFoldIdx) == m_QueryChainCount);
	}

void SCOP40Bench::StoreScore(uint ChainIndex1, uint ChainIndex2, float ScoreAB)
	{
	if (!KeepScore(ScoreAB))
		return;

	uint DomIdx1 = m_DomIdxs[ChainIndex1];
	uint DomIdx2 = m_DomIdxs[ChainIndex2];

	m_Scores.push_back(ScoreAB);
	m_DomIdx1s.push_back(DomIdx1);
	m_DomIdx2s.push_back(DomIdx2);
	}

float SCOP40Bench::AlignDomPair(uint ThreadIndex,
  uint Dom1, uint Dom2, uint &Lo1, uint &Lo2, string &Path)
	{
	Lo1 = UINT_MAX;
	Lo2 = UINT_MAX;
	Path.clear();

	const PDBChain &Chain1 = GetChainByDomIdx(Dom1);
	const PDBChain &Chain2 = GetChainByDomIdx(Dom2);

	const vector<vector<byte> > &Profile1 = GetProfileByDomIdx(Dom1);
	const vector<vector<byte> > &Profile2 = GetProfileByDomIdx(Dom2);

	asserta(ThreadIndex < SIZE(m_DAs));
	DSSAligner &DA = *m_DAs[ThreadIndex];
	DA.SetQuery(Chain1, &Profile1, 0, 0);
	DA.SetTarget(Chain2, &Profile2, 0, 0);
	DA.Align_NoAccel();
	Lo1 = DA.m_LoA;
	Lo2 = DA.m_LoB;
	Path = DA.m_PathAB;
	return DA.m_EvalueAB;
	}

void SCOP40Bench::OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA)
	{
	StoreScore(ChainIndex1, ChainIndex2, DA.m_EvalueAB);
	}

void SCOP40Bench::OnAlnBA(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA)
	{
	StoreScore(ChainIndex2, ChainIndex1, DA.m_EvalueBA);
	}

void SCOP40Bench::SetSFIdxToDomIdxs()
	{
	uint SFCount = GetSFCount();
	uint DomCount = GetDomCount();
	m_SFIdxToDomIdxs.clear();
	m_SFIdxToDomIdxs.resize(SFCount);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint SFIdx = m_DomIdxToSFIdx[DomIdx];
		m_SFIdxToDomIdxs[SFIdx].push_back(DomIdx);
		}
	}

void SCOP40Bench::SetDomIdxToL()
	{
	m_DomIdxToL.clear();
	uint DomCount = GetDomCount();
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		const string &Label = m_Doms[DomIdx];
		vector<string> Fields;
		Split(Label, Fields, '/');
		const string &Dom = Fields[0];
		map<string, uint>::const_iterator iter =
		  m_DomToChainIdx.find(Dom);
		asserta(iter != m_DomToChainIdx.end());
		const PDBChain &Chain = *m_Chains[iter->second];
		m_DomIdxToL.push_back(Chain.GetSeqLength());
		}
	}

void SCOP40Bench::SetDomIdxToHitIdxs()
	{
	uint DomCount = GetDomCount();
	uint HitCount = GetHitCount();
	m_DomIdxToHitIdxs.clear();
	m_DomIdxToHitIdxs.resize(DomCount);
	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		uint DomIdx1 = m_DomIdx1s[HitIdx];
		m_DomIdxToHitIdxs[DomIdx1].push_back(HitIdx);
		}
	}

void SCOP40Bench::ClearHits()
	{
	m_Scores.clear();
	m_DomIdx1s.clear();
	m_DomIdx2s.clear();
	}

float SCOP40Bench::GetVeryBadScore() const
	{
	if (m_ScoresAreEvalues)
		return 999.9f;
	else
		return -999.9f;
	}

bool SCOP40Bench::HitIsBetter(uint HitIdx1, uint HitIdx2) const
	{
	float Score1 = m_Scores[HitIdx1];
	float Score2 = m_Scores[HitIdx2];
	return ScoreIsBetter(Score1, Score2);
	}

bool SCOP40Bench::ScoreIsBetter(float Score1, float Score2) const
	{
	if (m_ScoresAreEvalues)
		return Score1 < Score2;
	else
		return Score1 > Score2;
	}

float SCOP40Bench::GetMeanLength(uint SFIdx) const
	{
	float Sum = 0;
	const vector<uint> &DomIdxs = m_SFIdxToDomIdxs[SFIdx];
	const uint N = SIZE(DomIdxs);
	asserta(N > 0);
	for (uint i = 0; i < N; ++i)
		{
		uint DomIdx = DomIdxs[i];
		const string &Label = m_Doms[DomIdx];
		vector<string> Fields;
		Split(Label, Fields, '/');
		const uint n = SIZE(Fields);
		asserta(n < 3);
		string Dom = Fields[0];
		map<string, uint>::const_iterator iter =
		  m_DomToChainIdx.find(Dom);
		asserta(iter != m_DomToChainIdx.end());
		uint ChainIdx = iter->second;
		const PDBChain &Chain = *m_Chains[ChainIdx];
		Sum += Chain.GetSeqLength();
		}
	return Sum/N;
	}

void SCOP40Bench::ScanDomHits()
	{
	m_ConsideredHitCount = 0;
	m_IgnoredHitCount = 0;
	m_DomIdxToScoreFirstFP.clear();
	m_DomIdxToSens1FP.clear();

	m_DomIdxToHitIdxFirstFP.clear();
	const uint DomCount = GetDomCount();
	const uint HitCount = GetHitCount();

	const float InitScore = GetVeryBadScore();
	m_DomIdxToScoreFirstFP.resize(DomCount, InitScore);
	m_DomIdxToSens1FP.resize(DomCount, 0);
	m_DomIdxToTP1Count.resize(DomCount, 0);
	m_DomIdxToHitIdxFirstFP.resize(DomCount, UINT_MAX);

	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		uint Dom1 = m_DomIdx1s[HitIdx];
		uint Dom2 = m_DomIdx2s[HitIdx];
		if (Dom1 == Dom2)
			continue;
		float Score = m_Scores[HitIdx];
		bool T = false;
		uint SF1 = m_DomIdxToSFIdx[Dom1];
		uint SF2 = m_DomIdxToSFIdx[Dom2];
		uint Fold1 = m_DomIdxToFoldIdx[Dom1];
		uint Fold2 = m_DomIdxToFoldIdx[Dom2];
		if (m_Mode == "family")
			T = (SF1 == SF2);
		else if (m_Mode == "fold")
			T = (Fold1 == Fold2);
		else if (m_Mode == "ignore")
			{
			if (Fold1 == Fold2 && SF1 != SF2)
				{
				++m_IgnoredHitCount;
				continue;
				}
			T = (Fold1 == Fold2);
			}
		else
			Die("m_Mode=%s", m_Mode.c_str());

		++m_ConsideredHitCount;
		if (T)
			{
			if (ScoreIsBetter(Score, m_DomIdxToScoreFirstFP[Dom1]))
				++m_DomIdxToTP1Count[Dom1];
			}
		else
			{
			if (ScoreIsBetter(Score, m_DomIdxToScoreFirstFP[Dom1]))
				{
				m_DomIdxToScoreFirstFP[Dom1] = Score;
				m_DomIdxToHitIdxFirstFP[Dom1] = HitIdx;
				}
			}
		}

	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		uint Dom1 = m_DomIdx1s[HitIdx];
		uint Dom2 = m_DomIdx2s[HitIdx];
		if (Dom1 == Dom2)
			continue;
		float Score = m_Scores[HitIdx];
		uint SF1 = m_DomIdxToSFIdx[Dom1];
		uint SF2 = m_DomIdxToSFIdx[Dom2];
		uint Fold1 = m_DomIdxToFoldIdx[Dom1];
		uint Fold2 = m_DomIdxToFoldIdx[Dom2];
		bool T = false;
		if (m_Mode == "family")
			T = (SF1 == SF2);
		else if (m_Mode == "fold")
			T = (Fold1 == Fold2);
		else if (m_Mode == "ignore")
			{
			if (Fold1 == Fold2 && SF1 != SF2)
				{
				++m_IgnoredHitCount;
				continue;
				}
			T = (Fold1 == Fold2);
			}
		else
			asserta(false);

		if (T)
			{
			if (ScoreIsBetter(Score, m_DomIdxToScoreFirstFP[Dom1]))
				m_DomIdxToSens1FP[Dom1] += 1;
			}
		}

	m_DomsWithHomologCount = 0;
	for (uint Dom = 0; Dom < DomCount; ++Dom)
		{
		uint SF = m_DomIdxToSFIdx[Dom];
		uint SFSize = m_SFSizes[SF];
		if (SFSize > 1)
			{
			++m_DomsWithHomologCount;
			if (m_DomIdxToTP1Count[Dom] > 0)
				++m_DomsWithHomologAndTP1Count;
			}
		}
	}

void SCOP40Bench::GetTPs1FP(vector<uint> &Doms1, vector<uint> &Doms2)
	{
	Doms1.clear();
	Doms2.clear();
	ScanDomHits();
	const uint HitCount = GetHitCount();

	uint GoodCount = 0;
	for (uint i = 0; i < HitCount; ++i)
		{
		uint Dom1 = m_DomIdx1s[i];
		uint Dom2 = m_DomIdx2s[i];
		float Score = m_Scores[i];
		uint SF1 = m_DomIdxToSFIdx[Dom1];
		uint SF2 = m_DomIdxToSFIdx[Dom2];
		if (Dom1 != Dom2 && SF1 == SF2)
			{
			Doms1.push_back(Dom1);
			Doms2.push_back(Dom2);
			}
		}
	}

uint SCOP40Bench::GetSens1stFP()
	{
	ScanDomHits();
	const uint HitCount = GetHitCount();

	uint GoodCount = 0;
	for (uint i = 0; i < HitCount; ++i)
		{
		uint Dom1 = m_DomIdx1s[i];
		uint Dom2 = m_DomIdx2s[i];
		if (Dom1 != Dom2 && IsT(Dom1, Dom2) == 1)
			{
			float Score = m_Scores[i];
			if (ScoreIsBetter(Score, m_DomIdxToScoreFirstFP[Dom1]))
				++GoodCount;
			}
		}
	return GoodCount;
	}

void SCOP40Bench::WriteBit(const string &FileName) const
	{
	if (FileName == "")
		return;
	const uint DomCount = SIZE(m_DomIdxToSFIdx);
	const uint SFCount = SIZE(m_SFs);
	const uint HitCount = SIZE(m_DomIdx1s);
	ProgressLog("Write %u doms (%u fams), %s hits to %s\n",
	  DomCount, SFCount, IntToStr(HitCount), FileName.c_str());
	asserta(SIZE(m_DomIdx2s) == HitCount);
	asserta(SIZE(m_Scores) == HitCount);

	FILE *f = CreateStdioFile(FileName);
	WriteStdioFile(f, &DomCount, sizeof(DomCount));
	WriteStdioFile(f, &HitCount, sizeof(HitCount));
	WriteStdioFile(f, m_DomIdxToSFIdx.data(), DomCount*sizeof(uint));
	WriteStdioFile(f, m_DomIdx1s.data(), HitCount*sizeof(uint));
	WriteStdioFile(f, m_DomIdx2s.data(), HitCount*sizeof(uint));
	WriteStdioFile(f, m_Scores.data(), HitCount*sizeof(float));
	CloseStdioFile(f);
	}

void SCOP40Bench::SetStats(float MaxFPR)
	{
	SetTFs();

	GetROCSteps(m_ROCStepScores, m_ROCStepNTPs, m_ROCStepNFPs);

	SmoothROCSteps(m_ROCStepScores, m_ROCStepNTPs, m_ROCStepNFPs, 100, MaxFPR,
	  m_SmoothScores, m_SmoothNTPs, m_SmoothNFPs, m_SmoothTPRs, m_SmoothFPRs);

	ROCStepsToTsv(opt_roc, m_SmoothScores, m_SmoothNTPs, m_SmoothNFPs,
	  m_SmoothTPRs, m_SmoothFPRs);

	m_nt_epq0_1 = GetNTPAtEPQThreshold(m_ROCStepNTPs, m_ROCStepNFPs, 0.1f);
	m_nt_epq1 = GetNTPAtEPQThreshold(m_ROCStepNTPs, m_ROCStepNFPs, 1);
	m_nt_epq10 = GetNTPAtEPQThreshold(m_ROCStepNTPs, m_ROCStepNFPs, 10);
	m_nt_firstfp = GetSens1stFP();
	}

void SCOP40Bench::ReadHits_Tsv_DSS()
	{
	const string TsvFN = "d:/int/dss/out/scop40_dss.tsv";
	FILE *f = OpenStdioFile(TsvFN);
	uint FileSize = GetStdioFileSize32(f);
	Progress("Reading hits DSS\n");
	uint HitCount = 0;
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		if (HitCount%500000 == 0)
			{
			uint Pos = GetStdioFilePos32(f);
			Progress("Hits %.2f%%  %s\r", GetPct(Pos, FileSize), IntToStr(HitCount));
			}
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		string Dom1;
		string Dom2;
		GetDomFromLabel(Fields[1], Dom1);
		GetDomFromLabel(Fields[2], Dom2);
		const map<string, uint>::iterator iter1 = m_DomToIdx.find(Dom1);
		const map<string, uint>::iterator iter2 = m_DomToIdx.find(Dom2);
		if (iter1 == m_DomToIdx.end() || iter2 == m_DomToIdx.end())
			continue;
		uint DomIdx1 = iter1->second;
		uint DomIdx2 = iter2->second;
		float Score = (float) StrToFloat(Fields[0]);
		m_DomIdx1s.push_back(DomIdx1);
		m_DomIdx2s.push_back(DomIdx2);
		m_Scores.push_back(Score);
		if (DomIdx1 != DomIdx2)
			{
			m_DomIdx1s.push_back(DomIdx2);
			m_DomIdx2s.push_back(DomIdx1);
			m_Scores.push_back(Score);
			}
		++HitCount;
		}
	Progress("%u hits DSS\n", HitCount);
	}

void SCOP40Bench::ReadHits_Tsv(const string &Algo)
	{
	if (Algo == "dss")
		{
		ReadHits_Tsv_DSS();
		return;
		}
	const string TsvFN = "c:/data/scop40pdb/alignResults/rawoutput/" + Algo + "aln";
	uint ScoreFieldNr = 2;
	if (Algo == "foldseek" || Algo == "blastp")
		ScoreFieldNr = 10;
	FILE *f = OpenStdioFile(TsvFN);
	uint FileSize = GetStdioFileSize32(f);
	Progress("Reading hits %s\n", Algo.c_str());
	uint HitCount = 0;
	string Line;
	vector<string> Fields;
	uint BadLineCount = 0;
	while (ReadLineStdioFile(f, Line))
		{
		if (HitCount > 0 && HitCount%500000 == 0)
			{
			uint Pos = GetStdioFilePos32(f);
			Progress("Hits %.2f%%  %s\r", GetPct(Pos, FileSize), IntToStr(HitCount));
			}
		for (uint i = 0; i < SIZE(Line); ++i)
			if (Line[i] == ' ')
				Line[i] = '\t';
		Split(Line, Fields, '\t');
	// 3dblastswaln has blank line
		uint FieldCount = SIZE(Fields);
		if (FieldCount <= ScoreFieldNr)
			{
			++BadLineCount;
			continue;
			}
		string Label1 = Fields[0];
		string Label2 = Fields[1];
		string Dom1, Dom2;
		GetDomFromLabel(Label1, Dom1);
		GetDomFromLabel(Label2, Dom2);

		const map<string, uint>::iterator iter1 = m_DomToIdx.find(Dom1);
		const map<string, uint>::iterator iter2 = m_DomToIdx.find(Dom2);
		if (iter1 == m_DomToIdx.end() || iter2 == m_DomToIdx.end())
			continue;
		uint DomIdx1 = iter1->second;
		uint DomIdx2 = iter2->second;
		asserta(FieldCount > ScoreFieldNr);
		float Score = (float) StrToFloat(Fields[ScoreFieldNr]);
		m_DomIdx1s.push_back(DomIdx1);
		m_DomIdx2s.push_back(DomIdx2);
		m_Scores.push_back(Score);
		++HitCount;
		}
	ProgressLog("%u hits, %u bad lines %s\n", HitCount, BadLineCount, Algo.c_str());
	}

void SCOP40Bench::WriteOutputFiles()
	{
	if (optset_scoredist)
		ScoreDist(opt_scoredist);

	if (optset_epqx)
		{
		ProgressLog("Writing %s\n", opt_epqx);
		FILE *fRocx = CreateStdioFile(opt_epqx);
		EPQToTsvX(fRocx, 100);
		CloseStdioFile(fRocx);
		}
	if (optset_tps1fp)
		{
		FILE *f = CreateStdioFile(opt_tps1fp);
		vector<uint> Doms1;
		vector<uint> Doms2;
		GetTPs1FP(Doms1, Doms2);
		const uint N = SIZE(Doms1);
		asserta(SIZE(Doms2) == N);
		DSSAligner &DA = *m_DAs[0];
		for (uint i = 0; i < N; ++i)
			{
			ProgressStep(i, N, "Re-aligning TPs");

			uint Dom1 = Doms1[i];
			uint Dom2 = Doms2[i];
			if (Dom1 > Dom2)
				continue;
			string Path;
			uint Lo1, Lo2;
			float Evalue = AlignDomPair(0, Dom1, Dom2,
			  Lo1, Lo2, Path);
			if (optset_evalue && Evalue > opt_evalue)
				continue;

			string CIGAR;
			LocalPathToCIGAR(Path.c_str(), Lo1, Lo2,  CIGAR, true);
			const char *DomName1 = m_Doms[Dom1].c_str();
			const char *DomName2 = m_Doms[Dom2].c_str();
			fprintf(f, "T\t%.3g\t%s\t%s\t%s\n",
			  Evalue, DomName1, DomName2, CIGAR.c_str());
			}
		if (optset_nfp)
			{
			ProgressLog("Realigning FPs...");
			vector<uint> FPHitIdxs;
			uint HitCount = GetHitCount();
			for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
				{
				uint Dom1 = m_DomIdx1s[HitIdx];
				uint Dom2 = m_DomIdx2s[HitIdx];
				if (IsT(Dom1, Dom2) == 0)
					FPHitIdxs.push_back(HitIdx);
				}
			std::random_device rd;
			std::mt19937 g(rd());

			shuffle(FPHitIdxs.begin(), FPHitIdxs.end(), g);
			const uint n = min(opt_nfp, SIZE(FPHitIdxs));
			for (uint j = 0; j < n; ++j)
				{
				uint k = FPHitIdxs[j];
				uint Dom1 = m_DomIdx1s[k];
				uint Dom2 = m_DomIdx2s[k];

				string Path;
				uint Lo1, Lo2;
				float Evalue = AlignDomPair(0, Dom1, Dom2, 
				  Lo1, Lo2, Path);

				string CIGAR;
				LocalPathToCIGAR(Path.c_str(), Lo1, Lo2, CIGAR, true);
				const char *DomName1 = m_Doms[Dom1].c_str();
				const char *DomName2 = m_Doms[Dom2].c_str();
				fprintf(f, "F\t%.3g\t%s\t%s\t%s\n",
					Evalue, DomName1, DomName2, CIGAR.c_str());
				}
			ProgressLog(" done.\n");
			}
		CloseStdioFile(f);
		}
	}

void cmd_scop40bench()
	{
	string CalFN;
	if (g_Arg1 == ".")
#ifdef _MSC_VER
		CalFN = "c:/src/reseek_scop40/reseek_db/scop40_family.cal";
#else
		CalFN = "/c/src/reseek_scop40/reseek_db/scop40_family.cal";
#endif
	else
		CalFN = g_Arg1;
	SCOP40Bench SB;
	asserta(optset_benchmode);
	DSSParams Params;
	SB.ReadChains(CalFN, "");

	Params.SetFromCmdLine();
	Params.m_DBSize = (float) SB.m_ChainCount;

	SB.Setup(Params);
	if (optset_scores_are_not_evalues)
		SB.m_ScoresAreEvalues = false;
	float MaxFPR = 0.005f;
	if (optset_maxfpr)
		MaxFPR = (float) opt_maxfpr;
	ResetTimers();
	SB.m_QuerySelf = true;
	SB.m_ScoresAreEvalues = true;
	if (opt_scores_are_not_evalues)
		SB.m_ScoresAreEvalues = false;
	SB.Run();
	SB.WriteBit(opt_savebit);
	SB.SetStats(MaxFPR);
	SB.WriteOutputFiles();

	uint AlnCount = DSSAligner::m_AlnCount;
	uint HitCount = SB.GetHitCount();
	uint UFilterCount = DSSAligner::m_UFilterCount;
	uint ComboFilterCount = DSSAligner::m_ComboFilterCount;
	uint Secs = SB.m_Secs;
	float AlnsPerThreadPerSec = SB.m_AlnsPerThreadPerSec;
	uint nt_firstfp = SB.m_nt_firstfp;
	float SensFirstFP = float(nt_firstfp)/SB.m_NT;
	float SensEPQ0_1 = float(SB.m_nt_epq0_1)/SB.m_NT;
	float SensEPQ1 = float(SB.m_nt_epq1)/SB.m_NT;
	float SensEPQ10 = float(SB.m_nt_epq10)/SB.m_NT;
	uint TotalFilterCount = DSSAligner::m_UFilterCount + DSSAligner::m_ComboFilterCount;
	double FilterPct = GetPct(TotalFilterCount, AlnCount);
	double UFilterPct = GetPct(DSSAligner::m_UFilterCount, AlnCount);
	uint ComboFilterInputCount = AlnCount - DSSAligner::m_UFilterCount;
	double ComboFilterPct =
	  GetPct(DSSAligner::m_ComboFilterCount, ComboFilterInputCount);
	double FoundFract1 = double(SB.m_DomsWithHomologAndTP1Count)/
	  SB.m_DomsWithHomologCount;

	ProgressLog("SEPQ0.1=%.4f", SensEPQ0_1);
	ProgressLog(" SEPQ1=%.4f", SensEPQ1);
	ProgressLog(" SEPQ10=%.4f", SensEPQ10);
	ProgressLog(" S1FP=%.4f", SensFirstFP);
	ProgressLog(" FF1=%.4f", FoundFract1);
	ProgressLog(" secs=%u", Secs);
	ProgressLog(" mode=%s\n", SB.m_Mode.c_str());

	if (DSSAligner::m_UFilterCount > 0)
		Log("UFil=%.1f\n", UFilterPct);
	if (DSSAligner::m_ComboFilterCount > 0)
		Log(" CFil=%.1f\n", ComboFilterPct);
	if (DSSAligner::m_ParasailSaturateCount > 0)
		Log(" Sat=%u\n", DSSAligner::m_ParasailSaturateCount);

	Log("QP cache hits %u, misses %u\n",
	  SB.m_QPCacheHits.load(), SB.m_QPCacheMisses.load());

	Log("NT %u NF %u NI %u considered %u ignored %u\n",
	  SB.m_NT, SB.m_NF, SB.m_NI, SB.m_ConsideredHitCount, SB.m_IgnoredHitCount);
	}

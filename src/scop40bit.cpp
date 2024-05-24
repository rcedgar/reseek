#include "myutils.h"
#include "scop40bit.h"
#include "scop40bench.h"
#include "pdbchain.h"
#include "sort.h"

/***
==> 3dblastaln <==
d12asa_ d12asa_ 1981
d12asa_ d1nnha_ 428
d12asa_ d1eova2 406

==> 3dblastswaln <==
                    <<< blank line!!  
d12asa_ d12asa_ 1981
d12asa_ d1nnha_ 428

==> cealn <==
d12asa_ d12asa_ 8.03 0.00 1.00
d12asa_ d1nnha_ 6.81 2.40 0.81
d12asa_ d1b8aa2 6.70 2.76 0.75

==> cleswaln <==
d12asa_ d12asa_ 15284
d12asa_ d1b8aa2 5739
d12asa_ d1nnha_ 5620

==> dalialn <==
d12asa_ d12asa_ 58.0 0.0
d12asa_ d1nnha_ 27.0 2.7
d12asa_ d1b8aa2 25.1 2.7

==> foldseekaln <==
      0       1     2     3   4 5   6     7 8     9        10     11
d1a1xa_	d1a1xa_	0.000	106	105	0	1	106	1	106	2.549E-163	1045
d1a1xa_	d1jsga_	0.000	108	105	0	1	106	4	111	9.137E-93	622
d1a1xa_	d3saoa_	0.000	77	76	0	26	102	27	103	2.448E-21	188

==> mmseqsaln <==
d12asa_ d12asa_ 674
d12asa_ d2gz4a1 39
d12asa_ d1b8aa2 30

==> tmaln <==
d12asa_ d12asa_ 1.0
d12asa_ d1b8aa2 0.75576
d12asa_ d1nnha_ 0.751037

==> tmfastaln <==
d12asa_ d12asa_ 1.0000 1.0000
d12asa_ d1b8aa2 0.7558 0.7394
d12asa_ d1nnha_ 0.7510 0.8309
***/

const uint MAGIC1 = MAGIC('h', 'i', 't', '$');
const uint MAGIC2 = MAGIC('$', 'h', 'i', 't');

const uint MaxFamLen = 16;

void GetDomFamFoldFromDomSlashFam(const string &Label,
 string &Dom, string &Fam, string &Fold)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	if (SIZE(Fields) == 2)
		{
		Dom = Fields[0];
		Fam = Fields[1];
		Split(Fam, Fields, '.');
		asserta(SIZE(Fields) == 4);
		Fold = Fields[0] + '.' + Fields[1] + '.' + Fields[2];
		}
	else
		{
		Dom = Label;
		Fam = "-";
		Fold = "-";
		}
	}

void GetDomFamFromDomSlashFam(const string &Label,
 string &Dom, string &Fam)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	if (SIZE(Fields) == 2)
		{
		Dom = Fields[0];
		Fam = Fields[1];
		}
	else
		{
		Dom = Label;
		Fam = "-";
		}
	}

void GetDomFromDomSlashFam(const string &Label, string &Dom)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	//asserta(SIZE(Fields) == 2);
	Dom = Fields[0];
	}

void SCOP40Bit::ReadDomInfo(const string &ChainsFN)
	{
	string Line;
	vector<string> Fields;

	const string DomFamFN = "c:/int/scop40/out/domain_family.tsv";
	const string DomainsFN = "c:/int/scop40/out/algo_domains.txt";

	FILE *f1 = OpenStdioFile(DomainsFN);
	string Dom;
	while (ReadLineStdioFile(f1, Dom))
		{
		asserta(strlen(Dom.c_str()) == 7);
		uint Idx = SIZE(m_Doms);
		asserta(m_DomToIdx.find(Dom) == m_DomToIdx.end());
		m_Doms.push_back(Dom);
		m_DomToIdx[Dom] = Idx;
		}
	CloseStdioFile(f1);

	FILE *f2 = OpenStdioFile(DomFamFN);
	m_DomIdxToFamIdx.clear();
	const uint DomCount = SIZE(m_Doms);
	m_DomIdxToFamIdx.resize(DomCount, UINT_MAX);
	while (ReadLineStdioFile(f2, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		const string &Dom = Fields[0];
		const string &Fam = Fields[1];
		asserta(strlen(Fam.c_str()) < MaxFamLen);
		asserta(strlen(Dom.c_str()) == 7);
		map<string, uint>::iterator iter1 = m_DomToIdx.find(Dom);
		if (iter1 == m_DomToIdx.end())
			continue;
		uint DomIdx = iter1->second;
		uint FamIdx = UINT_MAX;
		map<string, uint>::const_iterator iter2 = m_FamToIdx.find(Fam);
		if (iter2 == m_FamToIdx.end())
			{
			FamIdx = SIZE(m_Fams);
			m_Fams.push_back(Fam);
			m_FamToIdx[Fam] = FamIdx;
			}
		else
			FamIdx = iter2->second;
		asserta(m_DomToFam.find(Dom) == m_DomToFam.end());
		m_DomToFam[Dom] = Fam;
		m_DomIdxToFamIdx[DomIdx] = FamIdx;
		}
	const uint FamCount = SIZE(m_Fams);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		asserta(m_DomIdxToFamIdx[DomIdx] < FamCount);
	CloseStdioFile(f2);
	}

static void ReadScop40Bit(const string &FN)
	{
	FILE *f = OpenStdioFile(FN);
	uint Bytes = GetStdioFileSize32(f);
	byte *Data = myalloc(byte, Bytes);
	ReadStdioFile(f, Data, Bytes);
	CloseStdioFile(f);

	const uint *uptr = (const uint *) Data;
	uint MAGIC = *uptr++;
	asserta(MAGIC == 0xD03FA3);
	uint DomCount = *uptr++;
	uint FamCount = *uptr++;
	ProgressLog("%u doms, %u fams, mfl %u\n",
	  DomCount, FamCount, MaxFamLen);
	vector<string> m_Doms;
	vector<string> m_Fams;
	vector<uint> DomIdxToFamIdx;
	const byte *bptr = (const byte *) uptr;
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		char Dom[8];
		for (uint i = 0; i < 8; ++i)
			Dom[i] = *bptr++;
		asserta(Dom[7] == 0);
		m_Doms.push_back(string(Dom));
		uptr = (const uint *) bptr;
		uint FamIdx = *uptr++;
		bptr = (const byte *) uptr;
		DomIdxToFamIdx.push_back(FamIdx);
		}
	asserta(SIZE(DomIdxToFamIdx) == DomCount);
	bptr = (const byte *) uptr;
	char *Fam = myalloc(char, MaxFamLen+1);
	for (uint FamIdx = 0; FamIdx < FamCount; ++FamIdx)
		{
		for (uint i = 0; i < MaxFamLen; ++i)
			{
			char c = *bptr++;
			Fam[i] = c;
			}
		Fam[MaxFamLen] = 0;
		string sFam = string(Fam);
		StripWhiteSpace(sFam);
		m_Fams.push_back(string(sFam));
		}
	uptr = (const uint *) bptr;
	uint HitCount = 0;
	while ((const byte *) uptr < Data + Bytes)
		{
		uint Idx1 = *uptr++;
		uint Idx2 = *uptr++;
		const float *fptr = (const float *) uptr;
		float Score = *fptr++;
		uptr = (const uint *) fptr;
		asserta(Idx1 < DomCount);
		asserta(Idx2 < DomCount);
		++HitCount;
		}
	ProgressLog("%s hits\n", IntToStr(HitCount));
	}

void SCOP40Bit::WriteHits_Bin(const string &FileName) const
	{
	if (FileName == "")
		return;
	const uint DomCount = SIZE(m_DomIdxToFamIdx);
	const uint FamCount = SIZE(m_Fams);
	const uint HitCount = SIZE(m_DomIdx1s);
	ProgressLog("Write %u doms (%u fams), %s hits to %s\n",
	  DomCount, FamCount, IntToStr(HitCount), FileName.c_str());
	asserta(SIZE(m_DomIdx2s) == HitCount);
	asserta(SIZE(m_Scores) == HitCount);

	FILE *f = CreateStdioFile(FileName);
	WriteStdioFile(f, &DomCount, sizeof(DomCount));
	WriteStdioFile(f, &HitCount, sizeof(HitCount));
	WriteStdioFile(f, m_DomIdxToFamIdx.data(), DomCount*sizeof(uint));
	WriteStdioFile(f, m_DomIdx1s.data(), HitCount*sizeof(uint));
	WriteStdioFile(f, m_DomIdx2s.data(), HitCount*sizeof(uint));
	WriteStdioFile(f, m_Scores.data(), HitCount*sizeof(float));
	CloseStdioFile(f);
	}

void SCOP40Bit::ReadHits_Bin(const string &FileName)
	{
	uint DomCount = UINT_MAX;
	uint HitCount = UINT_MAX;
	FILE *f = OpenStdioFile(FileName);
	ReadStdioFile(f, &DomCount, sizeof(DomCount));
	ReadStdioFile(f, &HitCount, sizeof(HitCount));
	Progress("%s hits, %u doms %s\n",
	  IntToStr(HitCount), DomCount, FileName.c_str());
	uint32 FileSize = GetStdioFileSize32(f);
	m_DomIdxToFamIdx.resize(DomCount, UINT_MAX);
	m_DomIdx1s.resize(HitCount, UINT_MAX);
	m_DomIdx2s.resize(HitCount, UINT_MAX);
	m_Scores.resize(HitCount, FLT_MAX);
	ReadStdioFile(f, m_DomIdxToFamIdx.data(), DomCount*sizeof(uint));
	ReadStdioFile(f, m_DomIdx1s.data(), HitCount*sizeof(uint));
	ReadStdioFile(f, m_DomIdx2s.data(), HitCount*sizeof(uint));
	ReadStdioFile(f, m_Scores.data(), HitCount*sizeof(float));
	CloseStdioFile(f);
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
		GetDomFromDomSlashFam(Fields[1], Dom1);
		GetDomFromDomSlashFam(Fields[2], Dom2);
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


void SCOP40Bit::ReadHits_Tsv_DSS()
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
		GetDomFromDomSlashFam(Fields[1], Dom1);
		GetDomFromDomSlashFam(Fields[2], Dom2);
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
		GetDomFromDomSlashFam(Label1, Dom1);
		GetDomFromDomSlashFam(Label2, Dom2);

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

void SCOP40Bit::ReadHits_Tsv(const string &Algo)
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
		GetDomFromDomSlashFam(Label1, Dom1);
		GetDomFromDomSlashFam(Label1, Dom2);

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

uint SCOP40Bit::GetDomCount() const
	{
	return SIZE(m_DomIdxToFamIdx);
	}

uint SCOP40Bit::GetFamCount() const
	{
	return SIZE(m_Fams);
	}

const char *SCOP40Bit::GetDomName(uint DomIdx) const
	{
	asserta(DomIdx < SIZE(m_Doms));
	return m_Doms[DomIdx].c_str();
	}

uint SCOP40Bit::GetDomIdx(const string &DomName) const
	{
	map<string, uint>::const_iterator iter =
	  m_DomToIdx.find(DomName);
	if (iter == m_DomToIdx.end())
		Die("GetDomIdx(%s)", DomName.c_str());
	return iter->second;
	}

uint SCOP40Bit::GetHitCount() const
	{
	uint HitCount = SIZE(m_DomIdx1s);
	asserta(SIZE(m_DomIdx2s) == HitCount);
	asserta(SIZE(m_Scores) == HitCount);
	return HitCount;
	}

void SCOP40Bit::GetFamSizes_Present(vector<uint> &FamSizes) const
	{
	FamSizes.clear();
	uint FamCount = GetFamCount();
	uint DomCount = GetDomCount();
	uint HitCount = GetHitCount();
	vector<bool> DomPresent(DomCount);
	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		uint DomIdx1 = m_DomIdx1s[HitIdx];
		uint DomIdx2 = m_DomIdx2s[HitIdx];
		DomPresent[DomIdx1] = true;
		DomPresent[DomIdx2] = true;
		}

	FamSizes.resize(FamCount, 0);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		if (DomPresent[DomIdx])
			{
			uint FamIdx = m_DomIdxToFamIdx[DomIdx];
			assert(FamIdx < FamCount);
			++FamSizes[FamIdx];
			}
		}
	}

void SCOP40Bit::GetFamSizes(vector<uint> &FamSizes) const
	{
	FamSizes.clear();
	uint FamCount = GetFamCount();
	uint DomCount = GetDomCount();
	FamSizes.resize(FamCount, 0);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint FamIdx = m_DomIdxToFamIdx[DomIdx];
		assert(FamIdx < FamCount);
		++FamSizes[FamIdx];
		}
	}

void SCOP40Bit::CalcNXs(uint &NT, uint &NF) const
	{
	NT = 0;
	NF = 0;
	const uint DomCount = GetDomCount();
	const uint FamCount = GetFamCount();

// order matters, include self hit
//   (query,ref) is not equivalent to (ref,query)
	uint TotalPairCount = DomCount*DomCount;

	vector<uint> FamSizes;
	GetFamSizes(FamSizes);
	for (uint FamIdx = 0; FamIdx < FamCount; ++FamIdx)
		{
		uint Size = FamSizes[FamIdx];
		uint FamPairCount = Size*Size;
		NT += FamPairCount;
		}
	NF = TotalPairCount - NT;
	}

void SCOP40Bit::CalcNXs_Present(uint &NT, uint &NF) const
	{
	NT = 0;
	NF = 0;
	const uint DomCount = GetDomCount();
	const uint FamCount = GetFamCount();

// order matters, include self hit
//   (query,ref) is not equivalent to (ref,query)
	uint TotalPairCount = DomCount*DomCount;

	vector<uint> FamSizes;
	GetFamSizes_Present(FamSizes);
	for (uint FamIdx = 0; FamIdx < FamCount; ++FamIdx)
		{
		uint Size = FamSizes[FamIdx];
		uint FamPairCount = Size*Size;
		NT += FamPairCount;
		}
	NF = TotalPairCount - NT;
	}

bool SCOP40Bit::IsT(uint DomIdx1, uint DomIdx2) const
	{
	assert(DomIdx1 < SIZE(m_DomIdxToFamIdx));
	assert(DomIdx2 < SIZE(m_DomIdxToFamIdx));
	uint FamIdx1 = m_DomIdxToFamIdx[DomIdx1];
	uint FamIdx2 = m_DomIdxToFamIdx[DomIdx2];
	return FamIdx1 == FamIdx2;
	}

void SCOP40Bit::GetTFs(vector<bool> &TFs) const
	{
	uint HitCount = GetHitCount();
	TFs.clear();
	TFs.reserve(HitCount);
	uint DomCount = SIZE(m_Doms);
	uint FamCount = SIZE(m_Fams);
	for (uint i = 0; i < HitCount; ++i)
		{
		uint DomIdx1 = m_DomIdx1s[i];
		uint DomIdx2 = m_DomIdx2s[i];
		double TF = IsT(DomIdx1, DomIdx2);
		TFs.push_back(TF);
		}
	}

void SCOP40Bit::GetOrder(vector<uint> &Order) const
	{
	const uint HitCount = GetHitCount();
	Order.resize(HitCount);
	if (m_ScoresAreEvalues)
		QuickSortOrder(m_Scores.data(), HitCount, Order.data());
	else
		QuickSortOrderDesc(m_Scores.data(), HitCount, Order.data());
	}

void SCOP40Bit::GetROCSteps(vector<double> &ScoreSteps,
   vector<uint> &TPCounts, vector<uint> &FPCounts) const
	{
	const uint SMOOTHN = 100;
	ScoreSteps.clear();
	TPCounts.clear();
	FPCounts.clear();
	const uint HitCount = GetHitCount();

	ProgressLogPrefix("GetTFs()\n");
	vector<bool> TFs;
	GetTFs(TFs);
	asserta(SIZE(TFs) == HitCount);

	ProgressLogPrefix("GetOrder()\n");
	vector<uint> Order;
	GetOrder(Order);
	asserta(SIZE(Order) == HitCount);

	vector<double> TmpScoreSteps;
	vector<uint> TmpTPCounts;
	vector<uint> TmpFPCounts;
	double CurrentScore = m_Scores[Order[0]];
	uint NTP = 0;
	uint NFP = 0;
	uint LastNTP = 0;
	ProgressLogPrefix("TmpSteps\n");
	for (uint k = 0; k < HitCount; ++k)
		{
		uint i = Order[k];
		double Score = m_Scores[i];
		if (Score != CurrentScore)
			{
			TmpScoreSteps.push_back(CurrentScore);
			TmpTPCounts.push_back(NTP);
			TmpFPCounts.push_back(NFP);
			CurrentScore = Score;
			}
		if (TFs[i])
			++NTP;
		else
			++NFP;
		}
	TmpScoreSteps.push_back(CurrentScore);
	TmpTPCounts.push_back(NTP);
	TmpFPCounts.push_back(NFP);
	const uint n = SIZE(TmpScoreSteps);
	asserta(SIZE(TmpTPCounts) == n);
	asserta(SIZE(TmpFPCounts) == n);
	if (n <= 2*SMOOTHN)
		{
		ScoreSteps = TmpScoreSteps;
		TPCounts = TmpTPCounts;
		FPCounts = TmpFPCounts;
		return;
		}
	ProgressLog("SmoothSteps\n");
	for (uint Bin = 0; Bin < SMOOTHN; ++Bin)
		{
		uint Idx = UINT_MAX;
		if (Bin == 0)
			Idx = 0;
		else if (Bin + 1 == SMOOTHN)
			Idx = n - 1;
		else
			{
			Idx = (Bin*n)/SMOOTHN;
			asserta(Idx > 0 && Idx < n - 1);
			}
		ScoreSteps.push_back(TmpScoreSteps[Idx]);
		TPCounts.push_back(TmpTPCounts[Idx]);
		FPCounts.push_back(TmpFPCounts[Idx]);
		}
	}

void cmd_scop40bit()
	{
	const string &Algo = g_Arg1;
	asserta(optset_input);
	asserta(optset_output);

	SCOP40Bench SB;
	asserta(optset_benchmode);
	SB.m_Mode = string(opt_benchmode);
	SB.ReadChains(opt_input);
	SB.BuildDomFamIndexesFromQueryChainLabels();
	SB.ReadHits_Tsv(Algo);
	SB.WriteBit(opt_output);
	}

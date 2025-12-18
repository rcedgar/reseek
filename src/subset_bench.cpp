#include "myutils.h"
#include "subset_bench.h"
#include "alpha.h"
#include "sort.h"

static const uint MAGIC1 = 0xd06e1;
static const uint MAGIC2 = 0xd06e2;
static const uint MAGIC3 = 0xd06e3;
static const uint MAGIC4 = 0xd06e4;

void SubPattern(const string &Pattern, const string &x,	string &s)
	{
	s.clear();
	bool Found = false;
	for (uint i = 0; i < SIZE(Pattern); ++i)
		{
		char c = Pattern[i];
		if (c == '@')
			{
			asserta(!Found);
			Found = true;
			s += x;
			}
		else
			s += c;
		}
	asserta(Found);
	}

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values)
	{
	Names.clear();
	Values.clear();

	vector<string> Fields;
	Split(VarStr, Fields, ';');

	const uint n = SIZE(Fields);
	for (uint i = 0; i < n; ++i)
		{
		const string &NameEqValue = Fields[i];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			Die("SubsetBench::ParseVarStr(%s) not name=value '%s'",
				VarStr.c_str(), Fields[i].c_str());
		const string &Name = Fields2[0];
		const string &ValueStr = Fields2[1];
		float Weight = StrToFloatf(ValueStr);
		Names.push_back(Name);
		Values.push_back(Weight);
		}
	}

void SubsetBench::AddDom(const string &Dom, const string &ScopId)
	{
	vector<string> Fields;
	Split(ScopId, Fields, '.');
	asserta(SIZE(Fields) == 4);
	const string SF = Fields[0] + "." + Fields[1] + "." + Fields[2];
	uint SFIdx = UINT_MAX;
	if (m_SFToIdx.find(SF) == m_SFToIdx.end())
		{
		SFIdx = SIZE(m_SFs);
		m_SFs.push_back(SF);
		m_SFToIdx[SF] = SFIdx;
		}
	else
		SFIdx = m_SFToIdx[SF];
	m_SFIdxs.push_back(SFIdx);
	asserta(SFIdx < SIZE(m_SFIdxToSize));
	m_SFIdxToSize[SFIdx] += 1;

	uint DomIdx = UINT_MAX;
	if (m_DomToIdx.find(Dom) != m_DomToIdx.end())
		{
		Die("Duplicate dom >%s", Dom.c_str());
		DomIdx = m_DomToIdx[Dom];
		}
	DomIdx = SIZE(m_Doms);
	m_Doms.push_back(Dom + "/" + SF);
	m_DomToIdx[Dom] = DomIdx;
	}

void SubsetBench::ReadLookup(const string &FileName)
	{
	m_SFIdxToSize.clear();
	m_SFIdxToSize.resize(2000, 0);

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	vector<string> Fields2;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		const string &Dom = Fields[0];
		const string &ScopId = Fields[1];
		AddDom(Dom, ScopId);
		}
	CloseStdioFile(f);

	m_NT = 0;
	const uint SFCount = SIZE(m_SFs);
	for (uint SFIdx = 0; SFIdx < SFCount; ++SFIdx)
		{
		uint Size = m_SFIdxToSize[SFIdx];
		asserta(Size > 0);
		m_NT += (Size*(Size - 1))/2;
		}
	m_NT *= 2;
	}

void SubsetBench::AllocDope(uint DopeSize)
	{
	asserta(m_DopeSize == 0);
	m_DomIdxQs = myalloc(uint16_t, DopeSize);
	m_DomIdxTs = myalloc(uint16_t, DopeSize);
	m_TPs = myalloc(bool, DopeSize);
	m_DopeSize = DopeSize;
	}

void SubsetBench::AllocHits()
	{
	asserta(m_Scores == 0);
	asserta(m_ScoreOrder == 0);
	m_Scores = myalloc(float, m_DopeSize);
	m_ScoreOrder = myalloc(uint, m_DopeSize);
	}

uint SubsetBench::GetDomIdx(const string &Label, bool ErrOk) const
	{
	string Dom;
	SCOP40Bench::GetDomFromLabel(Label, Dom);
	map<string, uint>::const_iterator iter = m_DomToIdx.find(Dom);
	if (iter == m_DomToIdx.end())
		{
		if (ErrOk)
			return UINT_MAX;
		Die("SubsetBench::GetDomIdx(%s)", Label.c_str());
		}
	uint Idx = iter->second;
	return Idx;
	}

uint SubsetBench::GetSFIdx(uint DomIdx) const
	{
	asserta(DomIdx < SIZE(m_SFIdxs));
	return m_SFIdxs[DomIdx];
	}

// 1=Q, 2=T, 3=Evalue
// Keep lower triangle only
void SubsetBench::MakeDopeFromHits(const string &FN)
	{
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	vector<string> LabelQs;
	vector<string> LabelTs;
	uint HitCount = 0;
	while (ReadLineStdioFile(f, Line))
		{
		if (++HitCount%100000 == 0)
			Progress("Hits %u\r", HitCount);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) >= 3);
		const string &LabelQ = Fields[0];
		const string &LabelT = Fields[1];
		if (LabelQ == LabelT)
			continue;

		string DomQ, DomT;
		SCOP40Bench::GetDomFromLabel(LabelQ, DomQ);
		SCOP40Bench::GetDomFromLabel(LabelT, DomT);

		uint DomIdxQ = GetDomIdx(DomQ, true);
		uint DomIdxT = GetDomIdx(DomT, true);
		if (DomIdxQ == UINT_MAX || DomIdxT == UINT_MAX)
			continue;

		double Evalue = StrToFloat(Fields[2]);
		if (Evalue >= 10)
			continue;
		if (LabelQ > LabelT)
			continue;

		LabelQs.push_back(LabelQ);
		LabelTs.push_back(LabelT);
		}
	Progress("Hits %u\n", HitCount);

	uint DopeSize = SIZE(LabelQs);
	asserta(SIZE(LabelTs) == DopeSize);
	ProgressLog("Dope size %u\n", DopeSize);
	AllocDope(DopeSize);

	uint NT = 0;
	uint NF = 0;
	for (uint Idx = 0; Idx < m_DopeSize; ++Idx)
		{
		const string &LabelQ = LabelQs[Idx];
		const string &LabelT = LabelTs[Idx];

		string DomQ, DomT;
		SCOP40Bench::GetDomFromLabel(LabelQ, DomQ);
		SCOP40Bench::GetDomFromLabel(LabelT, DomT);

		uint DomIdxQ = GetDomIdx(DomQ, false);
		uint DomIdxT = GetDomIdx(DomT, false);

		m_DopeDomIdxs.insert(DomIdxQ);
		m_DopeDomIdxs.insert(DomIdxT);

		m_DomIdxQs[Idx] = DomIdxQ;
		m_DomIdxTs[Idx] = DomIdxT;

		uint SFIdxQ = GetSFIdx(DomIdxQ);
		uint SFIdxT = GetSFIdx(DomIdxT);

		if (SFIdxQ == SFIdxT)
			{
			++NT;
			m_TPs[Idx] = true;
			}
		else
			{
			++NF;
			m_TPs[Idx] = false;
			}
		}

	CloseStdioFile(f);

	ProgressLog("%u / %u doms in dope, %u TPs, %u FPs\n",
		SIZE(m_DopeDomIdxs), SIZE(m_Doms), NT, NF);
	}

void SubsetBench::WriteDope(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);

	WriteStdioFile(f, &MAGIC1, sizeof(MAGIC1));
	WriteStdioFile(f, &m_DopeSize, sizeof(m_DopeSize));
	WriteStdioFile(f, (void *) m_DomIdxQs, m_DopeSize*sizeof(m_DomIdxQs[0]));
	WriteStdioFile(f, (void *) m_DomIdxTs, m_DopeSize*sizeof(m_DomIdxTs[0]));
	WriteStdioFile(f, (void *) m_TPs, m_DopeSize*sizeof(m_TPs[0]));
	WriteStdioFile(f, &MAGIC2, sizeof(MAGIC2));

	CloseStdioFile(f);
	}

void SubsetBench::ReadDope(const string &FN)
	{
	if (FN == "")
		return;
	FILE *f = OpenStdioFile(FN);

	uint Word;
	ReadStdioFile(f, &Word, sizeof(Word));
	asserta(Word == MAGIC1);

	uint DopeSize;
	ReadStdioFile(f, &DopeSize, sizeof(DopeSize));
	AllocDope(DopeSize);

	ReadStdioFile(f, (void *) m_DomIdxQs, m_DopeSize*sizeof(m_DomIdxQs[0]));
	ReadStdioFile(f, (void *) m_DomIdxTs, m_DopeSize*sizeof(m_DomIdxTs[0]));
	ReadStdioFile(f, (void *) m_TPs, m_DopeSize*sizeof(m_TPs[0]));

	ReadStdioFile(f, &Word, sizeof(Word));
	asserta(Word == MAGIC2);

	CloseStdioFile(f);
	}

void SubsetBench::MakeByteSeqs(
	const string &FN, BYTE_SEQ_FN BSFn, uint AlphaSize,
	const string &BSFN) const
	{
	if (BSFN == "")
		return;

	FILE *f = CreateStdioFile(BSFN);
	vector<vector<byte> > ByteSeqVec;

	const uint DomCount = SIZE(m_Doms);
	ByteSeqVec.clear();
	ByteSeqVec.resize(DomCount);

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	const uint ChainCount = SIZE(Chains);
	uint SumLength = 0;
	for (uint i = 0; i < ChainCount; ++i)
		{
		PDBChain &Chain = *Chains[i];
		const string &Label = Chain.m_Label;
		string Dom;
		SCOP40Bench::GetDomFromLabel(Label, Dom);
		map<string, uint>::const_iterator iter = m_DomToIdx.find(Dom);
		if (iter == m_DomToIdx.end())
			continue;
		uint DomIdx = iter->second;
		asserta(DomIdx < DomCount);
		vector<byte> &ByteSeq = ByteSeqVec[DomIdx];
		asserta(ByteSeq.empty());
		BSFn(Chain, AlphaSize, ByteSeq);
		uint L = SIZE(ByteSeq);
		asserta(L == Chain.GetSeqLength());
		SumLength += L;
		}

	for (uint i = 0; i < m_DopeSize; ++i)
		{
		uint DomIdxQ = m_DomIdxQs[i];
		uint DomIdxT = m_DomIdxTs[i];
		if (ByteSeqVec[DomIdxQ].empty())
			Die("Missing chain >%s", m_Doms[DomIdxQ].c_str());
		if (ByteSeqVec[DomIdxT].empty())
			Die("Missing chain >%s", m_Doms[DomIdxT].c_str());
		}

	byte *ByteSeqs = myalloc(byte, SumLength);
	uint Offset = 0;
	vector<uint> Ls;
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		const vector<byte> &ByteSeq = ByteSeqVec[DomIdx];
		const uint L = SIZE(ByteSeq);
		Ls.push_back(L);
		memcpy(ByteSeqs + Offset, ByteSeq.data(), L);
		Offset += L;
		}
	asserta(Offset == SumLength);

	WriteStdioFile(f, &MAGIC3, sizeof(MAGIC3));
	WriteStdioFile(f, &DomCount, sizeof(DomCount));
	WriteStdioFile(f, &SumLength, sizeof(SumLength));
	WriteStdioFile(f, (void *) Ls.data(), DomCount*sizeof(Ls[0]));
	WriteStdioFile(f, (void *) ByteSeqs, SumLength);
	WriteStdioFile(f, &MAGIC4, sizeof(MAGIC4));

	CloseStdioFile(f);
	}

void SubsetBench::ReadByteSeqs(const string &FN, uint AS,
	vector<vector<byte> > &ByteSeqs) const
	{
	FILE *f = OpenStdioFile(FN);

	uint Word;
	ReadStdioFile(f, &Word, sizeof(Word));
	asserta(Word == MAGIC3);

	uint DomCount, SumLength;
	ReadStdioFile(f, &DomCount, sizeof(DomCount));
	ReadStdioFile(f, &SumLength, sizeof(SumLength));

	asserta(DomCount == SIZE(m_Doms));
	uint *Ls = myalloc(uint, DomCount);
	byte *BSData = myalloc(byte, SumLength);

	ReadStdioFile(f, (void *) Ls, DomCount*sizeof(Ls[0]));
	ReadStdioFile(f, (void *) BSData, SumLength);

	ReadStdioFile(f, &Word, sizeof(Word));
	asserta(Word == MAGIC4);

	ByteSeqs.clear();
	ByteSeqs.resize(DomCount);
	uint Offset = 0;
	for (uint i = 0; i < DomCount; ++i)
		{
		uint L = Ls[i];
		vector<byte> &ByteSeq = ByteSeqs[i];
		ByteSeq.resize(L);
		//memcpy(ByteSeq.data(), BSData + Offset, L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			byte Letter = BSData[Offset++];
			asserta(Letter < AS);
			ByteSeq[Pos] = Letter;
			}
		}
	asserta(Offset == SumLength);

	CloseStdioFile(f);
	}

void SubsetBench::ByteSeqsToFasta(vector<vector<byte> > &ByteSeqs,
	const string &FN) const
	{
	FILE *f = CreateStdioFile(FN);
	const uint SeqCount = SIZE(ByteSeqs);
	for (uint DomIdx = 0; DomIdx < SeqCount; ++DomIdx)
		{
		const vector<byte> &ByteSeq = ByteSeqs[DomIdx];
		const uint L = SIZE(ByteSeq);
		string Seq;
		Seq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			byte Letter = ByteSeq[Pos];
			if (Letter == 0xff)
				Seq += 'A'; // UNDEFINED_ZERO_OVERLOAD
			else
				Seq += g_LetterToCharAmino[Letter];
			}
		const string &Label = m_Doms[DomIdx];
		SeqToFasta(f, Label, Seq);
		}
	CloseStdioFile(f);
	}

float **SubsetBench::AllocSWMx() const
	{
	asserta(m_MaxL > 0);
	float **SWMx = myalloc(float *, m_MaxL);
	for (uint i = 0; i < m_MaxL; ++i)
		SWMx[i] = myalloc(float, m_MaxL);
	return SWMx;
	}

void SubsetBench::LoadByteSeqs(const vector<string> &FNs)
	{
	const uint FeatureCount = SIZE(m_AlphaSizes);
	asserta(SIZE(FNs) == FeatureCount);
	asserta(m_ByteSeqVec.empty());
	asserta(m_ByteSeqFNs.empty());
	m_ByteSeqFNs = FNs;
	m_ByteSeqVec.resize(FeatureCount);
	uint DomCount = GetDomCount();
	for (uint i = 0; i < FeatureCount; ++i)
		{
		if (i == 0)
			m_MaxL = 0;
		const uint AS = m_AlphaSizes[i];
		const string &FN = FNs[i];
		Progress("Load %s\n", FN.c_str());
		vector<vector<byte> > &ByteSeqs = m_ByteSeqVec[i];
		ReadByteSeqs(FN, AS, ByteSeqs);
		asserta(SIZE(ByteSeqs) == DomCount);
		if (i == 0)
			{
			m_Ls.clear();
			m_Ls.reserve(DomCount);
			for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
				{
				uint L = SIZE(ByteSeqs[DomIdx]);
				m_MaxL = max(L, m_MaxL);
				m_Ls.push_back(L);
				}
			}
		else
			{
			for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
				asserta(SIZE(ByteSeqs[DomIdx]) == m_Ls[DomIdx]);
			}
		}
	}

float * const * SubsetBench::MakeSWMx(uint ThreadIdx,
	uint DomIdxQ, uint DomIdxT)
	{
	asserta(ThreadIdx < m_ThreadCount);

	const uint DomCount = SIZE(m_Ls);
	asserta(DomIdxQ < DomCount);
	asserta(DomIdxT < DomCount);

	const uint LQ = m_Ls[DomIdxQ];
	const uint LT = m_Ls[DomIdxT];
	asserta(LQ <= m_MaxL);
	asserta(LT <= m_MaxL);

	float **SWMx = m_SWMxs[ThreadIdx];
	const uint FeatureCount = GetFeatureCount();
	asserta(FeatureCount > 0);

	const float *SubstMxPtr0 = m_SubstMxPtrs[0];
	const float w0 = m_Weights[0];
	const uint AS0 = m_AlphaSizes[0];
	const byte *Q0 = m_ByteSeqVec[0][DomIdxQ].data();
	const byte *T0 = m_ByteSeqVec[0][DomIdxT].data();
	for (uint PosQ = 0; PosQ < LQ; ++PosQ)
		{
		byte q = Q0[PosQ];
		assert(q < AS0);
		for (uint PosT = 0; PosT < LT; ++PosT)
			{
			byte t = T0[PosT];
			assert(t < AS0);
			SWMx[PosQ][PosT] = w0*SubstMxPtr0[AS0*q + t];
			}
		}

	for (uint FeatureIdx = 1; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		const float *SubstMxPtr = m_SubstMxPtrs[FeatureIdx];
		const float w = m_Weights[FeatureIdx];
		const uint AS = m_AlphaSizes[FeatureIdx];
		const byte *Q = m_ByteSeqVec[FeatureIdx][DomIdxQ].data();
		const byte *T = m_ByteSeqVec[FeatureIdx][DomIdxT].data();
		for (uint PosQ = 0; PosQ < LQ; ++PosQ)
			{
			byte q = Q[PosQ];
			assert(q < AS);
			for (uint PosT = 0; PosT < LT; ++PosT)
				{
				byte t = T[PosT];
				assert(t < AS);
				SWMx[PosQ][PosT] += w*SubstMxPtr[AS*q + t];
				}
			}
		}
	return SWMx;
	}

void SubsetBench::ThreadBody(uint ThreadIdx)
	{
	asserta(ThreadIdx < SIZE(m_Mems));
	asserta(ThreadIdx < SIZE(m_SWMxs));
	XDPMem &Mem = *m_Mems[ThreadIdx];
	for (;;)
		{
		uint PairIdx = m_NextPairIdx++;
		if (PairIdx >= m_DopeSize)
			return;
		uint DomIdxQ = m_DomIdxQs[PairIdx];
		uint DomIdxT = m_DomIdxTs[PairIdx];
		float * const *SWMx = MakeSWMx(ThreadIdx, DomIdxQ, DomIdxT);
		uint LQ = m_Ls[DomIdxQ];
		uint LT = m_Ls[DomIdxT];
		assert(LQ <= m_MaxL);
		assert(LT <= m_MaxL);
		float Score = m_AF(Mem, LQ, LT, m_Open, m_Ext, SWMx);
		asserta(!isnan(Score));
		asserta(!isinf(Score));
		m_Scores[PairIdx] = Score;
		}
	}

void SubsetBench::Search(ALIGN_FN AF)
	{
	m_AF = AF;
	m_ThreadCount = GetRequestedThreadCount();
	if (m_Mems.empty())
		{
		for (uint ThreadIdx = 0; ThreadIdx < m_ThreadCount; ++ThreadIdx)
			{
			m_Mems.push_back(new XDPMem);
			m_SWMxs.push_back(AllocSWMx());
			}
		}
	else
		{
		asserta(SIZE(m_Mems) == m_ThreadCount);
		asserta(SIZE(m_SWMxs) == m_ThreadCount);
		}

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBody, this, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}

void SubsetBench::StaticThreadBody(SubsetBench *SB, uint ThreadIdx)
	{
	SB->ThreadBody(ThreadIdx);
	}

float *SubsetBench::ReadScoreMx(const string &FN, uint &AlphaSize) const
	{
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	FEATURE F = StrToFeature(Fields[0].c_str());
	AlphaSize = StrToUint(Fields[1]);
	float *ScoreMxPtr = myalloc(float, AlphaSize*AlphaSize);
	for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == AlphaSize+1);
		asserta(StrToUint(Fields[0]) == Letter);
		for (uint Letter2 = 0; Letter2 < AlphaSize; ++Letter2)
			ScoreMxPtr[Letter*AlphaSize + Letter2] =
				StrToFloatf(Fields[Letter2+1]);
		}
	CloseStdioFile(f);
	return ScoreMxPtr;
	}

void SubsetBench::SetWeights(const vector<float> &Weights)
	{
	const uint FeatureCount = GetFeatureCount();
	asserta(SIZE(Weights) == FeatureCount);
	float Sum = 0;
	for (uint i = 0; i < FeatureCount; ++i)
		Sum += Weights[i];
	asserta(Sum > 0);
	m_Weights.clear();
	for (uint i = 0; i < FeatureCount; ++i)
		m_Weights.push_back(Weights[i]/Sum);
	}

void SubsetBench::LoadScoreMxs(vector<string> &FNs)
	{
	const uint N = SIZE(FNs);
	asserta(m_SubstMxPtrs.empty());
	asserta(m_AlphaSizes.empty());
	m_SubstMxPtrs.resize(N);
	for (uint i = 0; i < N; ++i)
		{
		const string &FN = FNs[i];
		Progress("Load %s\n", FN.c_str());
		uint AS;
		m_SubstMxPtrs[i] = ReadScoreMx(FN, AS);
		m_AlphaSizes.push_back(AS);
		}
	}

void cmd_subset_bench_dope()
	{
	const string &HitsFN = g_Arg1;
	const string &LookupFN = opt(lookup);

	SubsetBench SB;
	SB.ReadLookup(LookupFN);
	SB.MakeDopeFromHits(HitsFN);
	SB.WriteDope(opt(output));

	if (0)
		{
		SubsetBench SB2;
		SB2.ReadDope(opt(output));
		SB2.WriteDope("dope.tmp");
		}
	}

static void ByteSeq_AA(const PDBChain &Chain,
	uint AlphaSize,
	vector<byte> &ByteSeq)
	{
	asserta(AlphaSize == 20);
	const string &Seq = Chain.m_Seq;
	const uint L = SIZE(Seq);
	ByteSeq.clear();
	ByteSeq.reserve(L);
	for (uint i = 0; i < L; ++i)
		{
		byte Letter = g_CharToLetterAmino[Seq[i]];
		if (Letter >= 20)
			Letter = UNDEFINED_ZERO_OVERLOAD;
		ByteSeq.push_back(Letter);
		}
	}

void cmd_subset_bench_bsaa()
	{
	const string &DopeFN = g_Arg1;
	const string &LookupFN = opt(lookup);

	SubsetBench SB;
	SB.ReadLookup(LookupFN);
	SB.ReadDope(DopeFN);
	SB.MakeByteSeqs(opt(input), ByteSeq_AA, 20, opt(output));
	}

static FEATURE s_F;

static void ByteSeq_Feature(const PDBChain &Chain, uint AS,
	vector<byte> &ByteSeq)
	{
	DSS D;
	D.Init(Chain);
	const uint L = Chain.GetSeqLength();
	ByteSeq.clear();
	ByteSeq.reserve(L);
	for (uint i = 0; i < L; ++i)
		{
		uint Letter = D.GetFeature(s_F, i);
		if (Letter >= AS)
			Letter = UNDEFINED_ZERO_OVERLOAD;
		ByteSeq.push_back(Letter);
		}
	}

void cmd_subset_bench_bsfeature()
	{
	const string &DopeFN = g_Arg1;
	const string &LookupFN = opt(lookup);
	asserta(optset_feature);
	const char *FeatureName = opt(feature);
	s_F = StrToFeature(FeatureName);
	uint AS = DSSParams::GetAlphaSize(s_F);

	SubsetBench SB;
	SB.ReadLookup(LookupFN);
	SB.ReadDope(DopeFN);
	SB.MakeByteSeqs(opt(input), ByteSeq_Feature, AS, opt(output));
	}

void cmd_subset_bench_bs2fa()
	{
	asserta(optset_output);
	asserta(optset_alpha_size);
	const string &DopeFN = g_Arg1;
	const string &LookupFN = opt(lookup);
	const string &BSFN = opt(input);
	const uint AS = opt(alpha_size);

	SubsetBench SB;
	SB.ReadLookup(LookupFN);
	SB.ReadDope(DopeFN);

	vector<vector<byte> > ByteSeqs;
	SB.ReadByteSeqs(BSFN, AS, ByteSeqs);
	SB.ByteSeqsToFasta(ByteSeqs, opt(output));
	}

void SubsetBench::SetScoreOrder()
	{
	if (m_ScoreOrder == 0)
		m_ScoreOrder = myalloc(uint, m_DopeSize);
	QuickSortOrderDesc(m_Scores, m_DopeSize, m_ScoreOrder);
	}

void SubsetBench::Bench(const string &Msg)
	{
	SetScoreOrder();
	asserta(m_ScoreOrder != 0);
	asserta(m_NT > 0);
	uint K = m_DopeSize;
	uint nt = 0;
	uint nf = 0;
	float LastScore = FLT_MAX;
	float SEPQ0_1 = FLT_MAX;
	float SEPQ1 = FLT_MAX;
	float SEPQ10 = FLT_MAX;
	const uint DomCount = GetDomCount();
	for (uint k = 0; k < K; ++k)
		{
		uint HitIdx = m_ScoreOrder[k];
		float Score = m_Scores[HitIdx];
		if (Score != LastScore)
			{
			if (Score >= LastScore)
				Die("k=%u HitIdx=%u Score=%.3g LastScore=%3g\n",
					k, HitIdx, Score, LastScore);
			float EPQ = 2*float(nf)/DomCount;
			float Sens = 2*float(nt)/m_NT;
			if (SEPQ0_1 == FLT_MAX && EPQ >= 0.1) SEPQ0_1 = Sens;
			if (SEPQ1 == FLT_MAX   && EPQ >= 1)   SEPQ1   = Sens;
			if (SEPQ10 == FLT_MAX  && EPQ >= 10)  SEPQ10  = Sens;
			LastScore = Score;
			}
		if (m_TPs[HitIdx])
			++nt;
		else
			++nf;
		}
	float EPQ = 2*float(nf)/DomCount;
	float Sens = 2*float(nt)/m_NT;
	if (SEPQ0_1 == FLT_MAX) SEPQ0_1 = Sens;
	if (SEPQ1 == FLT_MAX)   SEPQ1   = Sens;
	if (SEPQ10 == FLT_MAX)  SEPQ10  = Sens;
	m_Sum3 = SEPQ0_1*2 + SEPQ1*3/2 + SEPQ10;

	if (Msg != "")
		ProgressLog("%s ", Msg.c_str());
	ProgressLog("SEPQ0.1=%.3f", SEPQ0_1);
	ProgressLog(" SEPQ1=%.3f", SEPQ1);
	ProgressLog(" SEPQ10=%.3f", SEPQ10);
	ProgressLog(" Sum3=%.3f", m_Sum3);
	ProgressLog("\n");
	}

static float AF(
	XDPMem &Mem,
	uint LQ, uint LT,
	float Open, float Ext,
	float * const * SWMx)
	{
	float SWFast(XDPMem &Mem, const float * const *SMxData, uint LA, uint LB,
	  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	  string &Path);

	uint Loi, Loj, Leni, Lenj;
	string Path;
	float Score = SWFast(Mem, SWMx, LQ, LT, -Open, -Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

void cmd_subset_bench()
	{
	asserta(optset_bspattern);
	asserta(optset_mxpattern);
	asserta(optset_varstr);
	const string &DopeFN = g_Arg1;
	const string &LookupFN = opt(lookup);
	const string &BSPattern = opt(bspattern);
	const string &MxPattern = opt(mxpattern);
	const string &VarStr = opt(varstr);

	vector<string> FeatureNames;
	vector<float> Weights;
	ParseVarStr(VarStr, FeatureNames, Weights);

	vector<string> BSFNs;
	vector<string> SMFNs;
	const uint n = SIZE(FeatureNames);
	asserta(SIZE(Weights) == n);
	for (uint i = 0; i < n; ++i)
		{
		const string &Name = FeatureNames[i];
		string BSFN, MxFN;
		SubPattern(BSPattern, Name, BSFN);
		SubPattern(MxPattern, Name, MxFN);
		BSFNs.push_back(BSFN);
		SMFNs.push_back(MxFN);
		}

	SubsetBench SB;
	SB.ReadLookup(LookupFN);
	SB.ReadDope(DopeFN);
	SB.AllocHits();

	//SMFNs.push_back("../2025-12-16_train_alphadata_h/AA_20.out");
	//BSFNs.push_back("scop40c.bsaaf");
	//Weights.push_back(1);

	SB.m_Open = 3;
	SB.m_Ext = 0.3f;
	SB.LoadScoreMxs(SMFNs);
	SB.LoadByteSeqs(BSFNs);
	SB.SetWeights(Weights);
	SB.Search(AF);
	SB.Bench();
	}

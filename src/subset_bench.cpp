#include "myutils.h"
#include "subset_bench.h"
#include "alpha.h"

static const uint MAGIC1 = 0xd06e1;
static const uint MAGIC2 = 0xd06e2;
static const uint MAGIC3 = 0xd06e3;
static const uint MAGIC4 = 0xd06e4;

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
	}

void SubsetBench::AllocDope(uint DopeSize)
	{
	asserta(m_DopeSize == 0);
	m_DomIdxQs = myalloc(uint16_t, DopeSize);
	m_DomIdxTs = myalloc(uint16_t, DopeSize);
	m_TPs = myalloc(bool, DopeSize);
	m_DopeSize = DopeSize;
	}

uint SubsetBench::GetDomIdx(const string &Label) const
	{
	string Dom;
	SCOP40Bench::GetDomFromLabel(Label, Dom);
	map<string, uint>::const_iterator iter = m_DomToIdx.find(Dom);
	asserta(iter != m_DomToIdx.end());
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

	for (uint Idx = 0; Idx < m_DopeSize; ++Idx)
		{
		const string &LabelQ = LabelQs[Idx];
		const string &LabelT = LabelTs[Idx];

		uint DomIdxQ = GetDomIdx(LabelQ);
		uint DomIdxT = GetDomIdx(LabelT);

		m_DopeDomIdxs.insert(DomIdxQ);
		m_DopeDomIdxs.insert(DomIdxT);

		m_DomIdxQs[Idx] = DomIdxQ;
		m_DomIdxTs[Idx] = DomIdxT;

		uint SFIdxQ = GetSFIdx(DomIdxQ);
		uint SFIdxT = GetSFIdx(DomIdxT);

		if (SFIdxQ == SFIdxT)
			m_TPs[Idx] = true;
		else
			m_TPs[Idx] = false;
		}

	CloseStdioFile(f);

	ProgressLog("%u / %u doms in dope\n",
		SIZE(m_DopeDomIdxs), SIZE(m_Doms));
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
	const string &FN, BYTE_SEQ_FN BSFn, const string &BSFN)
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
		BSFn(Chain, ByteSeq);
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
	vector<byte> &ByteSeq)
	{
	const string &Seq = Chain.m_Seq;
	const uint L = SIZE(Seq);
	ByteSeq.clear();
	ByteSeq.reserve(L);
	for (uint i = 0; i < L; ++i)
		{
		byte Letter = g_CharToLetterAmino[Seq[i]];
		if (Letter >= 20)
			Letter = 0xff;
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
	SB.MakeByteSeqs(opt(input), ByteSeq_AA, opt(output));
	}

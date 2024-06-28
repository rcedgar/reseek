#include "myutils.h"
#include "pdbchain.h"
#include "bcadata.h"

void BCAData::Close()
	{
	if (m_Reading && !m_Writing)
		CloseReader();
	else if (m_Writing && !m_Reading)
		CloseWriter();
	else
		Die("BCAData::Close(), not open");
	}

void BCAData::Create(const string &FN)
	{
	if (FN == "")
		Die("Empty BCA filename");
	asserta(!m_Writing && !m_Reading);
	m_FN = FN;
	m_f = CreateStdioFile(FN);
	WriteStdioFile(m_f, &BCA_MAGIC, sizeof(BCA_MAGIC));

// Placeholder #1 will be overwritten with number of chains
// Placeholder #2 will be overwritten with address of labels in Close()
// Placeholder #3 will be overwritten with size of labels data
	uint64_t Placeholder = 0;
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	WriteStdioFile(m_f, &Placeholder, sizeof(Placeholder));
	m_Writing = true;
	}

void BCAData::WriteChain(const PDBChain &Chain)
	{
	asserta(m_Writing && !m_Reading);
	uint64_t Offset = GetStdioFilePos64(m_f);
	size_t n = m_Offsets.size();
	asserta(m_SeqLengths.size() == n);
	uint L = Chain.GetSeqLength();
	if (n > 0)
		{
		uint Ln_1 = m_SeqLengths[n-1];
		asserta(Offset == m_Offsets[n-1] + 7*Ln_1);
		}
	const string &Seq = Chain.m_Seq;
	uint Idx = SIZE(m_Labels);
	asserta(SIZE(m_SeqLengths) == Idx);

	m_Labels.push_back(Chain.m_Label);
	m_SeqLengths.push_back(L);
	m_Offsets.push_back(Offset);
	vector<uint16_t> ICs;
	Chain.GetICs(ICs);
	asserta(SIZE(ICs) == 3*L);
	WriteStdioFile64(m_f, Seq.c_str(), L);
	WriteStdioFile64(m_f, ICs.data(), 6*L);
	}

void BCAData::Open(const string &FN)
	{
	if (FN == "")
		Die("Empty BCA filename");
	asserta(!m_Writing && !m_Reading);
	asserta(m_f == 0);

	m_f = OpenStdioFile(FN);

	uint32_t Magic;
	ReadStdioFile(m_f, &Magic, sizeof(Magic));
	if (Magic != BCA_MAGIC)
		Die("Bad magic %08lx, invalid .bca file '%s'",
		  Magic, FN.c_str());

// Placeholder #1 will be overwritten with number of chains
// Placeholder #2 will be overwritten with address of labels in Close()
// Placeholder #3 will be overwritten with size of labels data
	uint64_t ChainCount64;
	ReadStdioFile(m_f, &ChainCount64, sizeof(uint64_t));
	ReadStdioFile(m_f, &m_SeqLengthsPos64, sizeof(uint64_t));
	ReadStdioFile(m_f, &m_LabelDataSize64, sizeof(uint64_t));
	uint64 Offset = GetStdioFilePos64(m_f);

	uint ChainCount = uint(ChainCount64);
	asserta(ChainCount == ChainCount64);

	m_SeqLengths.resize(ChainCount64);

	SetStdioFilePos64(m_f, m_SeqLengthsPos64);
	ReadStdioFile64(m_f, m_SeqLengths.data(), sizeof(uint32_t)*ChainCount);
	for (uint64 i = 0; i < ChainCount64; ++i)
		{
		m_Offsets.push_back(Offset);
		Offset += 7*m_SeqLengths[i];
		}

	uint LabelDataSize = uint(m_LabelDataSize64);
	asserta(LabelDataSize == m_LabelDataSize64);
	char *LabelData = myalloc(char, LabelDataSize);
	ReadStdioFile(m_f, LabelData, LabelDataSize);
	m_Labels.clear();
	uint n = 0;
	m_Labels.push_back(LabelData);
	for (uint i = 0; i + 1 < LabelDataSize; ++i)
		if (LabelData[i] == 0)
			m_Labels.push_back(LabelData + i + 1);
	uint LabelCount = SIZE(m_Labels);
	if (LabelCount != ChainCount64)
		Die("Bad BCA file, %u chains %u labels",
		  ChainCount, LabelCount);

	m_Reading = true;
	}

void BCAData::CloseReader()
	{
	asserta(m_Reading && !m_Writing);
	CloseStdioFile(m_f);
	m_f = 0;
	m_Reading = false;
	m_Writing = false;
	}

void BCAData::CloseWriter()
	{
	asserta(m_Writing && !m_Reading);

	const uint ChainCount = GetChainCount();
	m_SeqLengthsPos64 = GetStdioFilePos64(m_f);
	WriteStdioFile64(m_f, m_SeqLengths.data(), sizeof(uint32_t)*ChainCount);

	m_LabelDataSize64 = 0;
	for (uint i = 0; i < ChainCount; ++i)
		{
		const string &Label = m_Labels[i];
		const uint n = SIZE(Label) + 1;
		WriteStdioFile(m_f, Label.c_str(), n);
		m_LabelDataSize64 += n;
		}

// Re-wind to overwrite placeholders in file header
	SetStdioFilePos64(m_f, sizeof(BCA_MAGIC));
	uint64 ChainCount64 = ChainCount;

// #1 number of chains
// #2 address of labels
// #3 size of labels data
	WriteStdioFile(m_f, &ChainCount64, sizeof(ChainCount64));
	WriteStdioFile(m_f, &m_SeqLengthsPos64, sizeof(m_SeqLengthsPos64));
	WriteStdioFile(m_f, &m_LabelDataSize64, sizeof(m_LabelDataSize64));

	CloseStdioFile(m_f);
	m_f = 0;
	m_Writing = false;
	m_Reading = false;
	}

uint64 BCAData::GetICsOffset(uint64 ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_Offsets));
	asserta(ChainIdx < SIZE(m_SeqLengths));
	uint64 Offset = m_Offsets[ChainIdx];
	uint L = m_SeqLengths[ChainIdx];
	return Offset + L;
	}

uint64 BCAData::GetSeqOffset(uint64 ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_Offsets));
	return m_Offsets[ChainIdx];
	}

uint BCAData::GetSeqLength(uint64 ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_SeqLengths));
	return m_SeqLengths[ChainIdx];
	}

void BCAData::ReadChain(uint64 ChainIdx, PDBChain &Chain) const
	{
	asserta(m_Reading && !m_Writing);
	Chain.Clear();
	uint64 SeqOffset = GetSeqOffset(ChainIdx);
	uint L = GetSeqLength(ChainIdx);
	char *Seq = myalloc(char, L+1);
	ReadStdioFile64(m_f, SeqOffset, Seq, L);
	Seq[L] = 0;
	Chain.m_Seq = string(Seq);
	myfree(Seq);

	uint16_t *ICs = myalloc(uint16_t, 3*L);
	ReadStdioFile64(m_f, SeqOffset + L, ICs, 6*L);
	Chain.CoordsFromICs(ICs, L);
	asserta(ChainIdx < SIZE(m_Labels));
	Chain.m_Label = m_Labels[ChainIdx];
	}

void cmd_bca_stats()
	{
	BCAData BCA;
	BCA.Open(g_Arg1);
	uint ChainCount = BCA.GetChainCount();
	ProgressLog("%10u  Chains\n", ChainCount);
	uint64 SumL = 0;
	for (uint i = 0; i < ChainCount; ++i)
		SumL += BCA.m_SeqLengths[i];
	ProgressLog("%10u  Residues\n", SumL);
	ProgressLog("%10u  Label data bytes\n", (uint) BCA.m_LabelDataSize64);
	}

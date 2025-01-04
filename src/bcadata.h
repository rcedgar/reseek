#pragma once

class PDBChain;

// Binary C-alpha
class BCAData
	{
public:
	vector<string> m_Labels;
	vector<uint64_t> m_Offsets; // start of IC vector in file
	vector<uint32_t> m_SeqLengths;
	string m_FN;
	FILE *m_f = 0;
	bool m_Writing = false;
	bool m_Reading = false;
	uint64 m_SeqLengthsPos64 = UINT64_MAX;
	uint64 m_LabelDataSize64 = UINT64_MAX;

public:
	void Clear();
	void Create(const string &FN);
	void Open(const string &FN);
	void WriteChain(const PDBChain &Chain);
	void ReadChain(uint64 ChainIdx, PDBChain &Chain) const;
	void Close();
	uint GetChainCount() const { return SIZE(m_Labels); }
	uint64 GetSeqOffset(uint64 ChainIdx) const;
	//uint64 GetICsOffset(uint64 ChainIdx) const;
	uint GetSeqLength(uint64 ChainIdx) const;

private:
	void CloseWriter();
	void CloseReader();
	};

const uint32_t BCA_MAGIC = 0xBCABCA;

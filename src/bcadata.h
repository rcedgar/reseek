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

public:
	void Create(const string &FN);
	void Open(const string &FN);
	void WriteChain(const PDBChain &Chain);
	uint GetChainCount() const { return SIZE(m_Labels); }
	void Close();

private:
	void CloseWriter();
	void CloseReader();
	};

const uint32_t BCA_MAGIC = 0xBCABCA;

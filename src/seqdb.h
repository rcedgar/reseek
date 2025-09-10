#ifndef seqdb_h
#define seqdb_h

#include <map>
#include <set>

class SeqDB
	{
public:
	bool m_IsAligned;
	bool m_IsNucleo;
	bool m_IsNucleoSet;
	unsigned m_ColCount;
	vector<string> m_Labels;
	vector<string> m_Seqs;
	map<string, uint> m_LabelToIndex;

public:
	SeqDB()
		{
		m_IsAligned = false;
		m_IsNucleo = false;
		m_IsNucleoSet = false;
		m_ColCount = UINT_MAX;
		}

	void SetLabelToIndex();

	void ToLetters(const byte *CharToLetter);
	bool GetSeqByLabel(const string &Label, string &Seq,
	  bool FailOnError = true) const;
	uint GetSeqIndex(const string &Label, bool FailOnError = true) const;
	unsigned AddSeq(const string &Label, const string &Seq);
	const string &GetSeq(unsigned SeqIndex) const;
	const byte *GetByteSeq(uint SeqIndex) const;
	void GetSeq_StripGaps(unsigned SeqIndex, string &Seq, bool ToUpper = false) const;
	const string &GetLabel(unsigned SeqIndex) const;
	unsigned GetSeqLength(unsigned SeqIndex) const;
	uint GetUpperCount(unsigned uColIndex) const;
	uint GetLowerCount(unsigned uColIndex) const;
	uint GetLetterCount(unsigned uColIndex) const;
	uint GetGapCount(unsigned uColIndex) const;

	bool IsAligned() const;
	unsigned GetColCount() const;
	bool GetIsNucleo();
	unsigned GetSeqCount() const { return SIZE(m_Seqs); }
	void FromFasta(const string &FileName, bool AllowGaps = false);
	void FromFasta_Seqs(const string &FileName, const SeqDB &Seqs, bool AllowGaps = false);
	void WritePretty(FILE *f) const;
	void LogMe() const;
	void TruncLabels();

private:
	void SetIsNucleo();
	};

#endif // seqdb_h

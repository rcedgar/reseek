#ifndef sfasta_h
#define sfasta_h

#include "myutils.h"
#include <map>

class SeqInfo;

// Sequential reader for FASTA file format.
// Serves sequences in file order to save memory.
// Caches biggish chunks to compromise memory vs. speed.
class SFasta
	{
public:
	string m_FileName;
	FILE *m_File;
	bool m_AllowGaps;

	uint64 m_FileSize;

// Position to start next read
	uint64 m_FilePos;

// Cached data.
	char *m_Buffer;

// Bytes allocated to m_Buffer
	uint32 m_BufferSize;

// Current position in buffer, normally points to '>'
	uint32 m_BufferOffset;

// File data in buffer <= m_BufferSize
	uint32 m_BufferBytes;

// Current label
// Points into m_Buffer, not a separate buffer.
	char *m_Label;

// Current sequence length
	unsigned m_SeqLength;

// Current seq index
	unsigned m_SeqIndex;

	unsigned m_LongestLength;
	unsigned m_TooLongCount;

protected:
	bool m_IsNucleoSet;
	bool m_IsNucleo;
	map<char, unsigned> m_BadCharToCount;

public:
	SFasta();
	~SFasta();

	void Clear();
	void Open(const string &FileName);
	void Rewind();
	bool GetIsNucleo();

// Get next sequence.
// Returns zero on end-of-file
	const char *GetNextSeq();

// Length of most recent sequence returned by GetNextSeq().
	unsigned GetSeqLength() const { return m_SeqLength; }

// Label of most recent sequence returned by GetNextSeq().
	const char *GetLabel() const { return m_Label; }

// Index of most recent sequence returned by GetNextSeq().
	unsigned GetSeqIndex() const { return m_SeqIndex; }

	bool GetNextSI(SeqInfo &SI);

	unsigned GetPctDoneX10() const;
	double GetPctDone() const;

	void LogMe() const;

protected:
	void FillCache();
	const char *GetNextSeqLo();
	};

#endif // sfasta_h

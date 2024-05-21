#ifndef seqinfo_h
#define seqinfo_h

#include "obj.h"

class SeqInfo : public Obj
	{
	friend class ObjMgr;

public:
	unsigned m_Index;
	const char *m_Label;
	const byte *m_Seq;

// Buffers are non-zero iff memory owned by this object.
	char *m_LabelBuffer;
	byte *m_SeqBuffer;

	unsigned m_L;
	unsigned m_LabelBytes;
	unsigned m_MaxL;
	unsigned m_MaxLabelBytes;
	bool m_RevComp;
	bool m_IsORF;
	SeqInfo *m_ORFNucSeq;
	int m_ORFFrame;
	unsigned m_ORFNucLo;
	unsigned m_ORFNucHi;
	unsigned m_ORFNucL;

protected:
	SeqInfo(ObjMgr *OM);

public:
	~SeqInfo();
	virtual unsigned GetMemBytes() const
		{
		return m_MaxLabelBytes + m_MaxL;
		}
	virtual void OnZeroRefCount();

	unsigned GetIL() const
		{
		if (m_IsORF)
			return m_ORFNucL;
		return m_L;
		}

	void ToFasta(FILE *f) const
		{
		void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
		SeqToFasta(f, m_Seq, m_L, m_Label);
		}

	void ToFasta(const string &FileName) const
		{
		FILE *f = CreateStdioFile(FileName);
		ToFasta(f);
		CloseStdioFile(f);
		}

	void LogMe() const;
	void Init(unsigned Index);
	void AllocLabel(unsigned n);
	void AllocL(unsigned n);
	void AppendSeq(byte c);
	void AppendLabel(char c);
	void SetCopy(unsigned Index, const char *Label, const byte *Seq, unsigned L);
	void SetPtrs(unsigned Index, const char *Label, const byte *Seq, unsigned L);
	void Copy(const SeqInfo &rhs);
	void SetLabel(const char *Label);
	void GetRevComp(SeqInfo *RC);
	};

#endif // seqinfo_h

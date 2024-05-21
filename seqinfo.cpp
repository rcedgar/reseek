#include "myutils.h"
#include "seqinfo.h"
#include "alpha.h"
#include "objmgr.h"

SeqInfo::SeqInfo(ObjMgr *OM)  : Obj(OT_SeqInfo, OM)
	{
	m_Index = UINT_MAX;
	m_Label = 0;
	m_Seq = 0;
	m_LabelBuffer = 0;
	m_SeqBuffer = 0;
	m_L = 0;
	m_RevComp = false;
	m_LabelBytes = 0;
	m_MaxL = 0;
	m_MaxLabelBytes = 0;
	m_ORFNucSeq = 0;
	m_IsORF = false;
	m_ORFNucLo = UINT_MAX;
	m_ORFNucHi = UINT_MAX;
	m_ORFNucL = UINT_MAX;
	m_ORFFrame = 0;
	}

SeqInfo::~SeqInfo()
	{
	if (m_SeqBuffer != 0)
		myfree(m_SeqBuffer);
	if (m_LabelBuffer != 0)
		myfree(m_LabelBuffer);
	}

void SeqInfo::OnZeroRefCount()
	{
	m_Index = UINT_MAX;
	m_Seq = 0;
	m_Label = 0;
	m_L = 0;
	m_RevComp = false;
	m_IsORF = false;
	}

void SeqInfo::Copy(const SeqInfo &rhs)
	{
	AllocLabel(rhs.m_LabelBytes);
	AllocL(rhs.m_L);

	m_Index = rhs.m_Index;
	m_L = rhs.m_L;
	m_LabelBytes = rhs.m_LabelBytes;

	memcpy(m_LabelBuffer, rhs.m_Label, rhs.m_LabelBytes);
	memcpy(m_SeqBuffer, rhs.m_Seq, rhs.m_L);

	m_Seq = m_SeqBuffer;
	m_Label = m_LabelBuffer;

	m_IsORF = rhs.m_IsORF;
	m_ORFFrame = rhs.m_ORFFrame;
	m_ORFNucLo = rhs.m_ORFNucLo;
	m_ORFNucHi = rhs.m_ORFNucHi;
	m_ORFNucL = rhs.m_ORFNucL;
	m_ORFNucSeq = rhs.m_ORFNucSeq;
//	m_ORFNucSeq->Up();
	m_ObjMgr->Up(m_ORFNucSeq);
	}

void SeqInfo::Init(unsigned Index)
	{
	m_Index = Index;
	m_L = 0;
	m_LabelBytes = 0;
	m_Seq = m_SeqBuffer;
	m_Label = m_LabelBuffer;
	m_IsORF = false;
	}

void SeqInfo::AllocLabel(unsigned n)
	{
	if (n <= m_MaxLabelBytes)
		return;

	unsigned NewMaxLabelBytes = n + 128;
	char *NewLabelBuffer = myalloc(char, NewMaxLabelBytes);
	myfree(m_LabelBuffer);
	m_LabelBuffer = NewLabelBuffer;
	m_Label = NewLabelBuffer;
	m_MaxLabelBytes = NewMaxLabelBytes;
	}

void SeqInfo::AllocL(unsigned n)
	{
	if (n < m_MaxL)
		{
		m_Seq = m_SeqBuffer;
		return;
		}
	
	unsigned NewMaxL = n + 1024;
	byte *NewSeqBuffer = myalloc(byte, NewMaxL);
	if (m_L > 0)
		memcpy(NewSeqBuffer, m_Seq, m_L);
	myfree(m_SeqBuffer);
	m_Seq = NewSeqBuffer;
	m_SeqBuffer = NewSeqBuffer;
	m_MaxL = NewMaxL;
	}

void SeqInfo::SetLabel(const char *Label)
	{
	unsigned n = (unsigned) strlen(Label) + 1;
	AllocLabel(n);
	m_LabelBytes = n;
	memcpy(m_LabelBuffer, Label, n);
	m_Label = m_LabelBuffer;
	}

void SeqInfo::AppendSeq(byte c)
	{
	AllocL(m_L + 1);
	m_SeqBuffer[m_L++] = c;
	}

void SeqInfo::AppendLabel(char c)
	{
	if (m_LabelBytes >= m_MaxLabelBytes)
		{
		unsigned NewMaxLabelBytes = m_LabelBytes + 128;
		char *NewLabelBuffer = myalloc(char, NewMaxLabelBytes);
		if (m_LabelBytes > 0)
			{
			asserta(m_Label != 0);
			memcpy(NewLabelBuffer, m_Label, m_LabelBytes);
			}
		myfree(m_LabelBuffer);
		m_LabelBuffer = NewLabelBuffer;
		m_Label = NewLabelBuffer;
		m_MaxLabelBytes = NewMaxLabelBytes;
		}

	m_LabelBuffer[m_LabelBytes++] = c;
	}

void SeqInfo::SetCopy(unsigned Index, const char *Label, const byte *Seq, unsigned L)
	{
	m_Index = Index;
	SetLabel(Label);
	AllocL(L);
	memcpy(m_SeqBuffer, Seq, L);
	m_Seq = m_SeqBuffer;
	m_L = L;
	m_IsORF = false;
	}

void SeqInfo::SetPtrs(unsigned Index, const char *Label, const byte *Seq, unsigned L)
	{
	m_Index = Index;
	m_Label = Label;
	m_Seq = Seq;
	m_L = L;
	m_IsORF = false;
	}

void SeqInfo::GetRevComp(SeqInfo *RCSI)
	{
	RCSI->AllocL(m_L);
	RCSI->SetLabel(m_Label);
	RCSI->m_Index = UINT_MAX;
	byte *RCSeq = RCSI->m_SeqBuffer;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		byte rc = g_CharToCompChar[c];
		if (rc == INVALID_CHAR)
			rc = c;
		RCSeq[m_L - i - 1] = rc;
		}
	RCSI->m_L = m_L;
	RCSI->m_RevComp = !m_RevComp;
	RCSI->m_Index = m_Index;
	}

void SeqInfo::LogMe() const
	{
	Log("SeqInfo(%p) Seq %lx, Buff %lx, L %u, MaxL %u >%s\n",
	  this,
	  m_Seq,
	  m_SeqBuffer,
	  m_L,
	  m_MaxL,
	  m_Label);
	Log("%*.*s\n", m_L, m_L, m_Seq);
	}

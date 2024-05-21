#include "myutils.h"
#include "mx.h"
#include "seqdb.h"
//#include "timing.h"
#include <mutex>

#define StartTimer(x)	/* empty */
#define EndTimer(x)	/* empty */

char ProbToChar(float p);

static mutex g_Lock;
//list<MxBase *> *MxBase::m_Matrices = 0;
//unsigned MxBase::m_AllocCount;
//unsigned MxBase::m_ZeroAllocCount;
//unsigned MxBase::m_GrowAllocCount;
//double MxBase::m_TotalBytes;
//double MxBase::m_MaxBytes;

static const char *LogizeStr(const char *s)
	{
	double d = atof(s);
	d = log(d);
	return TypeToStr<float>(float(d));
	}

static const char *ExpizeStr(const char *s)
	{
	double d = atof(s);
	d = exp(d);
	return TypeToStr<float>(float(d));
	}

void MxBase::Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
  const SeqDB *DB, unsigned IdA, unsigned IdB)
	{
	StartTimer(MxBase_Alloc);

	if (DB != 0)
		{
		asserta(IdA != UINT_MAX);
		asserta(IdB != UINT_MAX);
		asserta(RowCount >= DB->GetSeqLength(IdA) + 1);
		asserta(ColCount >= DB->GetSeqLength(IdB) + 1);
		}
	if (RowCount > m_AllocatedRowCount || ColCount > m_AllocatedColCount)
		{
		if (m_AllocatedRowCount > 0)
			{
			if (opt_logmemgrows)
				Log("MxBase::Alloc grow %s %u x %u -> %u x %u, %s bytes\n",
				  Name, m_AllocatedRowCount, m_AllocatedColCount,
				  RowCount, ColCount,
				  IntToStr(GetBytes()));
			}


		EndTimer(MxBase_Alloc);
		StartTimer(MxBase_FreeData);
		FreeData();
		EndTimer(MxBase_FreeData);
		StartTimer(MxBase_Alloc);

		unsigned N = max(RowCount + 16, m_AllocatedRowCount);
		unsigned M = max(ColCount + 16, m_AllocatedColCount);
//		N = max(N, M); // TODO: gets too big

		EndTimer(MxBase_Alloc);
		StartTimer(MxBase_AllocData);
		AllocData(N, M);
		EndTimer(MxBase_AllocData);
		StartTimer(MxBase_Alloc);
		}
	
	unsigned n = sizeof(m_Name)-1;
	strncpy(m_Name, Name, n);
	m_Name[n] = 0;
	m_RowCount = RowCount;
	m_ColCount = ColCount;
	m_SeqDB = DB;
	m_IdA = IdA;
	m_IdB = IdB;

	EndTimer(MxBase_Alloc);
	}

void MxBase::LogMe(bool WithData, int Opts) const
	{
	Log("\n");
	if (Opts & OPT_EXP)
		Log("Exp ");
	else if (Opts & OPT_LOG)
		Log("Log ");
	bool ZeroBased = ((Opts & OPT_ZERO_BASED) != 0);
	Log("%s(%p) Rows %u/%u, Cols %u/%u",
	  m_Name, this,
	  m_RowCount, m_AllocatedRowCount,
	  m_ColCount, m_AllocatedColCount);
	if (m_SeqDB != 0 && m_IdA != UINT_MAX)
		Log(", A=%s", m_SeqDB->GetLabel(m_IdA).c_str());
	Log("\n");
	if (!WithData || m_RowCount == 0 || m_ColCount == 0)
		return;

	const char *z = GetAsStr(0, 0);
	unsigned Width = (unsigned) strlen(z);
	unsigned Mod = 1;
	for (unsigned i = 0; i < Width; ++i)
		Mod *= 10;

	if (m_Alpha[0] != 0)
		{
		Log("// Alphabet=%s\n", m_Alpha);
		Log("//      ");
		unsigned n = (unsigned) strlen(m_Alpha);
		for (unsigned j = 0; j < n; ++j)
			Log(" %*c", Width, m_Alpha[j]);
		Log("\n");
		for (unsigned i = 0; i < n; ++i)
			{
			Log("/* %c */ {", m_Alpha[i]);
			unsigned ci = m_Alpha[i];
			for (unsigned j = 0; j < n; ++j)
				{
				unsigned cj = m_Alpha[j];
				Log("%s,", GetAsStr(ci, cj));
				}
			Log("},  // %c\n", m_Alpha[i]);
			}
		return;
		}
	else if (m_Alpha2[0] != 0)
		{
		unsigned n = (unsigned) strlen(m_Alpha2);
		Log("// Alphabet=%s\n", m_Alpha2);
		Log("//      ");
		for (unsigned j = 0; j < n; ++j)
			Log(" %*c", Width, m_Alpha2[j]);
		Log("\n");
		for (unsigned i = 0; i < n; ++i)
			{
			Log("/* %c */ {", m_Alpha2[i]);
			unsigned ci = m_Alpha2[i];
			for (unsigned j = 0; j < n; ++j)
				Log("%s,", GetAsStr(i, j));
			Log("},  // %c\n", m_Alpha2[i]);
			}
		return;
		}

	const byte *A = 0;
	const byte *B = 0;
	//if (m_SeqDB != 0 && m_IdA != UINT_MAX)
	//	A = m_SeqDB->GetSeq(m_IdA);

	if (B != 0)
		{
		if (A != 0)
			Log("  ");
		Log("%5.5s", "");
		if (ZeroBased)
			for (unsigned j = 0; j < m_ColCount; ++j)
				Log("%*c", Width, B[j]);
		else
			for (unsigned j = 0; j < m_ColCount; ++j)
				Log("%*c", Width, j == 0 ? ' ' : B[j-1]);
		Log("\n");
		}

	if (A != 0)
		Log("  ");
	Log("%5.5s", "");
	for (unsigned j = 0; j < m_ColCount; ++j)
		Log("%*u", Width, j%Mod);
	Log("\n");

	for (unsigned i = 0; i < m_RowCount; ++i)
		{
		if (A != 0)
			{
			if (ZeroBased)
				Log("%c ", A[i]);
			else
				Log("%c ", i == 0 ? ' ' : A[i-1]);
			}
		Log("%4u ", i);
		
		for (unsigned j = 0; j < m_ColCount; ++j)
			{
			const char *s = GetAsStr(i, j);
			if (Opts & OPT_LOG)
				s = LogizeStr(s);
			else if (Opts & OPT_EXP)
				s = ExpizeStr(s);
			Log("%s", s);
			}
		Log("\n");
		}
	}
static unsigned g_MatrixFileCount;

void MxBase::LogCounts()
	{
	}

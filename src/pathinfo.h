#ifndef pathinfo_h
#define pathinfo_h

#include "obj.h"

#define AllocId	ALLOCID_pathinfo_hdr

class PathInfo : public Obj
	{
	friend class ObjMgr;

public:
	unsigned m_BufferBytes;
	char *m_Buffer;
	char *m_Path;

protected:
	PathInfo(ObjMgr *OM) : Obj(OT_PathInfo, OM)
		{
		m_BufferBytes = 0;
		m_Buffer = 0;
		m_Path = 0;
		}

	virtual ~PathInfo()
		{
		myfree(m_Buffer);
		m_Buffer = 0;
		}

public:
	virtual void Reset() { m_Path = 0; }
	virtual unsigned GetMemBytes() const
		{
		return m_BufferBytes;
		}

public:
	void Alloc(unsigned LA, unsigned LB)
		{
		unsigned Bytes = LA + LB + 1;
		if (Bytes > m_BufferBytes)
			{
			myfree(m_Buffer);
			m_BufferBytes = Bytes + 128;
			m_Buffer = myalloc(char, m_BufferBytes);
			}
		m_Path = m_Buffer;
		}

	void SetEmpty()
		{
	// Prevent crash if never alloc'd
		if (m_Path == 0)
			Alloc(1024, 1024);
		m_Path[0] = 0;
		}

	const char *Get() const
		{
		return m_Path;
		}

	const char *GetPath() const
		{
		return m_Path;
		}

	char *GetFront()
		{
		asserta(m_BufferBytes > 0);
		return m_Buffer;
		}

	char *GetBack()
		{
		asserta(m_BufferBytes > 0);
		return m_Buffer + m_BufferBytes - 1;
		}

	unsigned GetColCount() const
		{
		if (m_Path == 0)
			return 0;
		return unsigned(strlen(m_Path));
		}

	unsigned GetCounts(unsigned &M, unsigned &D, unsigned &I) const
		{
		M = 0;
		D = 0;
		I = 0;
		asserta(m_Path != 0);
		for (const char *p = m_Path; *p; ++p)
			{
			char c = *p;
			if (c == 'M')
				++M;
			else if (c == 'D')
				++D;
			else if (c == 'I')
				++I;
			else
				asserta(false);
			}
		return M + D + I;
		}

	void Reverse()
		{
		unsigned L = GetColCount();
		for (unsigned i = 0; i < L/2; ++i)
			swap(m_Path[i], m_Path[L-i-1]);
		}
	};

#undef AllocId
#endif // pathinfo_h

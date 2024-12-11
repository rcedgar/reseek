#ifndef xdpmem_h
#define xdpmem_h

#include "mx.h"

static const float MINUS_INFINITY = -9e9f;
static const float UNINIT = -8e8f;
static const int MINUS_INFINITY_INT = INT_MIN;

typedef float fn_SubstScore(void *UserData, uint PosA, uint PosB);
typedef float (*ptr_fn_SubstScore)(void *UserData, uint PosA, uint PosB);

class ObjMgr;

class XDPMem
	{
public:
	unsigned m_LA;
	unsigned m_LB;
	Mx<byte> m_TBBit;
	byte *m_RevA;
	byte *m_RevB;
	float *m_Buffer1;
	float *m_Buffer2;
	int *m_Buffer1_Int;

public:
	XDPMem()
		{
		Clear(true);
		}

	~XDPMem()
		{
		Clear(false);
		}

	void Clear(bool ctor = false)
		{
		if (!ctor)
			{
			m_TBBit.Clear();
			myfree(m_Buffer1);
			myfree(m_Buffer2);
			myfree(m_Buffer1_Int);
			myfree(m_RevA);
			myfree(m_RevB);
			}
		
		m_RevA = 0;
		m_RevB = 0;
		m_Buffer1 = 0;
		m_Buffer2 = 0;
		m_Buffer1_Int = 0;
		}

	byte *GetRevA()
		{
		return m_RevA;
		}

	byte *GetRevB()
		{
		return m_RevB;
		}

	byte **GetTBBit()
		{
		return m_TBBit.GetData();
		}

	int *GetDPRow1Int()
		{
		return m_Buffer1_Int + 1;
		}

	float *GetDPRow1()
		{
		return m_Buffer1 + 1;
		}

	float *GetDPRow2()
		{
		return m_Buffer2 + 1;
		}

	void Alloc(unsigned LA, unsigned LB)
		{
		Clear();
		m_LA = LA;
		m_LB = LB;
		m_TBBit.Alloc(LA+8, LB+8);

		m_Buffer1 = myalloc(float, LB+8);
		m_Buffer1_Int = myalloc(int, m_LB+8);
		m_Buffer2 = myalloc(float, m_LB+8);
		m_RevA = myalloc(byte, LA+8);
		m_RevB = myalloc(byte, LB+8);
		}
	};

float SWFast(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);

float XDropFwd(XDPMem &Mem,
  float X, float Open, float Ext, 
  fn_SubstScore SubFn, void *UserData,
  uint LoA, uint aLA, uint LoB, uint aLB,
  string &Path);

float XDropBwd(XDPMem &Mem,
  float X, float Open, float Ext, 
  fn_SubstScore SubFn, void *UserData, 
  uint HiA, uint LA, uint HiB, uint LB,
  string &Path);

void MergeFwdBwd(uint LA, uint LB,
  uint FwdLoA, uint FwdLoB, const string &FwdPath,
  uint BwdHiA, uint BwdHiB, const string &BwdPath,
  uint &LoA, uint &LoB, uint &HiA, uint &HiB, string &Path);

#endif // xdpmem_h

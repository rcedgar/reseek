#include "myutils.h"
#include "dssaligner.h"

float ViterbiFastMem(XDPMem &Mem, uint LA, uint LB,
	fn_SubstScore SubFn, void *UserData, string &Path);

void DSSAligner::AlignQueryTarget_Global()
	{
	uint LA = m_ChainA->GetSeqLength();
	uint LB = m_ChainB->GetSeqLength();

	m_GlobalScore =
		ViterbiFastMem(m_Mem, LA, LB, 
					   StaticSubstScore, (void *) this, m_GlobalPath);

	//void LogAln(const char *A, const char *B, const char *Path, unsigned ColCount);
	//const char *A = m_ChainA->m_Seq.c_str();
	//const char *B = m_ChainB->m_Seq.c_str();
	//const char *Path = m_GlobalPath.c_str();
	//uint ColCount = SIZE(m_GlobalPath);
	//LogAln(A, B, Path, ColCount);
	m_LoA = 0;
	m_LoB = 0;
	m_Path = m_GlobalPath;
	m_AlnFwdScore = m_GlobalScore;
	}

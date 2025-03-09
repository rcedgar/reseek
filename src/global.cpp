#include "myutils.h"
#include "dssaligner.h"

float ViterbiFastMem(XDPMem &Mem, uint LA, uint LB,
	fn_SubstScore SubFn, void *UserData, string &Path);

void DSSAligner::AlignQueryTarget_Global()
	{
	ClearAlign();

	bool MuAccept = true;
	if (m_Params->m_Omega > 0)
		{
		++m_MuFilterInputCount;
		bool MuFilterOk = MuDPFilter();
		if (!MuFilterOk)
			{
			MuAccept = false;
			++m_MuFilterDiscardCount;
			return;
			}
		}

	uint LA = m_ChainA->GetSeqLength();
	uint LB = m_ChainB->GetSeqLength();

	m_GlobalScore =
		ViterbiFastMem(m_Mem, LA, LB, 
					   StaticSubstScore, (void *) this, m_GlobalPath);
	m_LoA = 0;
	m_LoB = 0;
	m_Path = m_GlobalPath;
	}

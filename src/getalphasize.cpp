#include "myutils.h"
#include "dssparams.h"

uint DSSParams::GetAlphaSize(FEATURE F, bool FailOk)
	{
	switch (F)
		{
	case FEATURE_AA:
		return 20;

	case FEATURE_SS:
	case FEATURE_NENSS:
	case FEATURE_RENSS:
	case FEATURE_NormDens4:
	case FEATURE_NENDist4:
	case FEATURE_RENDist4:
	case FEATURE_AA4:
		return 4;

	case FEATURE_SS3:
	case FEATURE_NENSS3:
	case FEATURE_RENSS3:
	case FEATURE_AA3:
		return 3;

	case FEATURE_Conf:
	case FEATURE_NENConf:
	case FEATURE_RENConf:
	case FEATURE_NormDens:
	case FEATURE_NENDist:
	case FEATURE_RENDist:
	case FEATURE_HelixDens:
	case FEATURE_StrandDens:
	case FEATURE_DstNxtHlx:
	case FEATURE_DstPrvHlx:
	case FEATURE_NX:
	case FEATURE_PMDist:
		return 16;

	case FEATURE_Mu:
		return 36;

	case FEATURE_ConfU:
		return 17;
		}
	if (!FailOk)
		Die("GetAlphaSize(%s)", FeatureToStr(F));
	return UINT_MAX;
	}

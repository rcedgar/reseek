#include "myutils.h"
#include "museqsource.h"
#include "pdbchain.h"
#include "dss.h"
#include "seqinfo.h"

static int FloatToInt(float x)
	{
	return int(2*x);
	}

static float IntToFloat(int x)
	{
	return float(x)/2;
	}

static int GetIntCAlphaDist()
	{
	return FloatToInt(3.81f);
	}

static void CCA(const PDBChain &Chain)
	{
	const uint L = Chain.GetSeqLength();
	float x0 = Chain.m_Xs[0];
	float y0 = Chain.m_Ys[0];
	float z0 = Chain.m_Xs[0];

	float est_x = x0;
	float est_y = y0;
	float est_z = z0;

	float prev_x = x0;
	float prev_y = y0;
	float prev_z = z0;

	for (uint i = 1; i < L; ++i)
		{
		float x = Chain.m_Xs[i];
		float y = Chain.m_Ys[i];
		float z = Chain.m_Zs[i];
	
		int ix = FloatToInt(x);
		int iy = FloatToInt(y);
		int iz = FloatToInt(z);
	
		int prev_ix = FloatToInt(prev_x);
		int prev_iy = FloatToInt(prev_y);
		int prev_iz = FloatToInt(prev_z);

		int diff_x = ix - prev_ix;
		int diff_y = iy - prev_iy;
		int diff_z = iz - prev_iz;

		float dx = x - prev_x;
		float dy = y - prev_y;
		float dz = z - prev_z;
		float d = sqrt(dx*dx + dy*dy + dz*dz);

		float est_dx = IntToFloat(ix - prev_ix);
		float est_dy = IntToFloat(iy - prev_iy);

		est_x += est_dx;
		est_y += est_dy;

	// sqrt(dx*dx + dy*dy + dz*dz) = 3.81
	// dz = sqrt(3.81^2 - dx*dx - dy*dy)
		float est_dz_squared = 3.81f*3.81f - est_dx*est_dx - est_dy*est_dy;
		if (est_dz_squared < 0)
			est_dz_squared = 0;
		float est_dz = sqrt(est_dz_squared);
		float est_z_plus = est_z + est_dz;
		float est_z_minus = est_z - est_dz;
		if (fabs(z - est_z_plus) < fabs(z - est_z_minus))
			est_z = est_z_plus;
		else
			est_z = est_z_minus;

		Log("x=%8.1f (%8.1f)", x, est_x);
		Log(" y=%8.1f (%8.1f)", y, est_y);
		Log(" z=%8.1f (%8.1f)", z, est_z);
		Log(" | ix=%5d", ix);
		Log("  iy=%5d", iy);
		Log("  (dx = %8.1f", dx);
		Log("  est_dx = %8.1f)", est_dx);
		Log(" dy = %8.1f", dy);
		Log(" est_dy = %8.1f", est_dy);
		Log(" d = %.1f", d);
		Log("\n");

		prev_x = x;
		prev_y = y;
		prev_z = z;
		}
	}

void cmd_cca()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		CCA(*Chain);
		}
	}

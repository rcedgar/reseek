#include "myutils.h"
#include "pdbchain.h"
// #include "outputfiles.h"
#include "abcxyz.h"
#include <map>

float g_DALI_D = 20.0f;
float g_DALI_d0 = 0.2f;
float g_DALI_Theta = 0.2f;

static const int TBLSZ = 100;
static double *WeightLookup;

/***
DaliLite v5
comparemodules.f, line 1436
===========================
		enveloperadius=20.0
		x=1/(enveloperadius*enveloperadius)
		do i=0,100
				wght(i)=exp(-x*i*i)
		end do
***/
static double Weight(double y)
	{
	int iy = int(y+0.5);
	if (iy < 0)
		iy = 0;
	if (iy >= TBLSZ)
		iy = TBLSZ-1;
	double w2 = WeightLookup[iy];
	return w2;
	}

static double Weight_NoLookup(double y)
	{
	//const double D = 20.0;
	const double x = 1.0 / (g_DALI_D * g_DALI_D);
	double w = exp(-x * y * y);
	return w;
	}

static bool InitWeightLookup()
	{
	WeightLookup = myalloc(double, TBLSZ);
	for (int i = 0; i < TBLSZ; ++i)
		{
		double y = double(i);
		double w = Weight_NoLookup(y);
		WeightLookup[i] = w;
		}
	return true;
	};
static bool InitWeightLookupDone = InitWeightLookup();

/***
comparemodules.f, line 1397
  a, b are integer distances in units of 1/10 Angstrom,
  so multiply by 10 to get Angstroms.
===========================
		function dpscorefun(a,b) result(s)
		implicit none
		include 'parsizes.for'
		real s
		integer*2 a,b
c
		real x,y,d0
		logical lela
		parameter(lela=.true.)
		parameter(d0=0.20)
c !!!   elastic uses weights !!!
		x=float(abs(a-b))/10
		if(lela) then
				y=float(a+b)/20
				if(y.gt.100) then
						s=0.0
				else
						if(y.gt.0) then
						  s=wght(nint(y))*(d0-x/y)
						else
						  s=wght(nint(y))*d0
						end if
				end if
		end if
***/
static double dpscorefun(double a, double b)
	{
	double Score = 0;
	double x = fabs(a - b);
	double y = (a + b) / 2;
	if (y > 100)
		Score = 0;
	else
		{
		if (y > 0)
			Score = Weight(y) * (g_DALI_d0 - x / y);
		else
			Score = Weight(y) * g_DALI_d0;
		}
	return Score;
	}

double GetDALIScore_OffDiag(const PDBChain &Q, const PDBChain &T,
	const vector<uint> &PosQs, const vector<uint> &PosTs)
	{
	const uint Lali = SIZE(PosQs);
	asserta(SIZE(PosTs) == Lali);

	double Sum = 0;
	for (uint i = 0; i < Lali; ++i)
		{
		uint PosQi = PosQs[i];
		uint PosTi = PosTs[i];
		for (uint j = 0; j < Lali; ++j)
			{
			if (i == j)
				continue;

			uint PosQj = PosQs[j];
			uint PosTj = PosTs[j];

			double dij_Q = Q.GetDist(PosQi, PosQj);
			double dij_T = T.GetDist(PosTi, PosTj);
			double x = dpscorefun(dij_Q, dij_T);
			Sum += x;
			}
		}
	return Sum;
	}

double GetDALIScore(const PDBChain &Q, const PDBChain &T,
	const vector<uint> &PosQs, const vector<uint> &PosTs)
	{
	const uint Lali = SIZE(PosQs);
	double OffDiag = GetDALIScore_OffDiag(Q, T, PosQs, PosTs);
	double Score = OffDiag + Lali*g_DALI_Theta;
	Score /= 100.0;
	return Score;
	}

double GetDALIScore_Path(const PDBChain &Q, const PDBChain &T,
	const string &Path, uint LoQ, uint LoT)
	{
	const uint Cols = SIZE(Path);
	vector<uint> PosQs;
	vector<uint> PosTs;
	uint PosQ = LoQ;
	uint PosT = LoT;
	for (uint Col = 0; Col < Cols; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			PosQs.push_back(PosQ);
			PosTs.push_back(PosT);
			++PosT;
			++PosQ;
			break;

		case 'D':
			++PosQ;
			break;

		case 'I':
			++PosT;
			break;
			}
		}
	uint QL = Q.GetSeqLength();
	uint TL = T.GetSeqLength();
	if (PosQ > QL)
		Die("GetDALIScore_Path() PosQ=%u QL=%u", PosQ, QL);
	if (PosT > TL)
		Die("GetDALIScore_Path() PosT=%u TL=%u", PosT, TL);
	return GetDALIScore(Q, T, PosQs, PosTs);
	}

double GetDALIScore_Path_OffDiag(const PDBChain &Q, const PDBChain &T,
	const string &Path, uint LoQ, uint LoT)
	{
	const uint Cols = SIZE(Path);
	vector<uint> PosQs;
	vector<uint> PosTs;
	uint PosQ = LoQ;
	uint PosT = LoT;
	for (uint Col = 0; Col < Cols; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
			case 'M':
			PosQs.push_back(PosQ);
			PosTs.push_back(PosT);
			++PosT;
			++PosQ;
			break;

			case 'D':
			++PosQ;
			break;

			case 'I':
			++PosT;
			break;
			}
		}
	uint QL = Q.GetSeqLength();
	uint TL = T.GetSeqLength();
	if (PosQ > QL)
		Die("GetDALIScore_Path_OffDiag() PosQ=%u QL=%u", PosQ, QL);
	if (PosT > TL)
		Die("GetDALIScore_Path_OffDiag() PosT=%u TL=%u", PosT, TL);
	return GetDALIScore_OffDiag(Q, T, PosQs, PosTs);
	}

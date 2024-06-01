#include "myutils.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include <map>
#include <set>
#include "timing.h"

struct XDATA
	{
	const vector<vector<byte> > *ptrProfile1 = 0;
	const vector<vector<byte> > *ptrProfile2 = 0;
	DSSAligner *ptrDSS = 0;
	};

static float SubFn(void *UserData, uint Pos1, uint Pos2)
	{
	XDATA *XD = (XDATA *) UserData;
	DSSAligner &D = *XD->ptrDSS;
	float Score = D.GetScorePosPair(*XD->ptrProfile1, *XD->ptrProfile2, Pos1, Pos2);
	return Score;
	}

void MakeKmerToCoords(const vector<uint> &Kmers,
  map<uint, vector<uint> > &KmerToCoords)
	{
	KmerToCoords.clear();
	const uint n = SIZE(Kmers);
	for (uint Pos = 0; Pos < n; ++Pos)
		{
		uint Kmer = Kmers[Pos];
		if (KmerToCoords.find(Kmer) != KmerToCoords.end())
			KmerToCoords[Kmer].push_back(Pos);
		else
			{
			vector<uint> v;
			v.push_back(Pos);
			KmerToCoords[Kmer] = v;
			}
		}
	}

uint DSSAligner::GetU(const vector<uint> &Kmers1, const vector<uint> &Kmers2) const
	{
	set<uint> Set1;
	for (uint i = 0; i < SIZE(Kmers1); ++i)
		{
		uint Kmer = Kmers1[i];
		if (Kmer != UINT_MAX)
			Set1.insert(Kmer);
		}
	uint U = 0;
	for (uint i = 0; i < SIZE(Kmers2); ++i)
		{
		uint Kmer = Kmers2[i];
		if (Kmer != UINT_MAX && Set1.find(Kmer) != Set1.end())
			++U;
		}
	return U;
	}

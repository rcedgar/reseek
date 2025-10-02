#include "myutils.h"
#include "dssaligner.h"

float GetSelfRevScore(DSSAligner &DA, const DSSParams &Params,
	const PDBChain &Chain,
	const vector<vector<byte> > &Profile,
	const vector<byte> *ptrMuLetters,
	const vector<uint> *ptrMuKmers)
	{
	if (opt(selfrev0))
		return 0;

	DSS D;
	PDBChain RevChain;
	Chain.GetReverse(RevChain);
	vector<vector<byte> > RevProfile;
	D.Init(RevChain, Params);
	D.GetProfile(RevProfile);

	DA.SetQuery(Chain, &Profile, ptrMuLetters, ptrMuKmers, FLT_MAX);

	DA.SetTarget(RevChain, &RevProfile, ptrMuLetters, ptrMuKmers, FLT_MAX);
	DA.AlignQueryTarget();
	return DA.m_AlnFwdScore;
	}

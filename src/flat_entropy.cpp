#include "myutils.h"
#include "flat_chain.h"

static float get_entropy(
	const vector<uint8_t> &codes, uint start, uint W)
	{
	map<uint8_t, uint> unique2count;
	for (uint pos = start; pos < start + W; ++pos)
		{
		uint8_t code = codes[pos];
		map<uint8_t, uint>::const_iterator iter =
			unique2count.find(code);
		if (iter == unique2count.end())
			unique2count[code] = 1;
		else
			unique2count[code] += 1;
		}

	float H = 0;
	for (map<uint8_t, uint>::const_iterator iter = unique2count.begin();
		iter != unique2count.end(); ++iter)
		{
		uint count = iter->second;
		float P = float(count)/W;
		H += -P*log(P);
		}
	return H;
	}

float flat_get_entropy(
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const uint8_t *profQ,
	const uint8_t *profT,
	uint nfeat, uint fi)
	{
	const uint pathlen = uint(path.size());

	vector<uint32_t> posQs;
	vector<uint32_t> posTs;

	void path2posvecs(
		const string &path,
		uint loQ, uint LQ,
		uint loT, uint LT,
		vector<uint> &posQs,
		vector<uint> &posTs);
	path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
	const uint ncol = uint(posQs.size());

	vector<uint8_t> codeQs;
	vector<uint8_t> codeTs;

	codeQs.reserve(ncol);
	codeTs.reserve(ncol);

	for (uint col = 0; col < ncol; ++col)
		{
		uint posQ = posQs[col];
		uint posT = posTs[col];
		uint8_t codeQ = profQ[fi*LQ + posQ];
		uint8_t codeT = profT[fi*LT + posT];

		codeQs.push_back(codeQ);
		codeTs.push_back(codeT);
		}

	const uint W = 10;//@@TODO param
	float HQ = 0;
	float HT = 0;
	for (uint start = 0; start + W <= ncol; ++start)
		{
		HQ += get_entropy(codeQs, start, W);
		HT += get_entropy(codeTs, start, W);
		}
	return HQ + HT;
	}

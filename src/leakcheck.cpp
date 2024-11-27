#include "myutils.h"
#include <unordered_map>
#include "sort.h"

#if LEAK_CHECK
extern unordered_map<void *, pair<string, size_t> > *g_AllocMap;
void LeakCheck(const char *Msg)
	{
	vector<string> Locvec;
	vector<uint> Bytesvec;
	for (unordered_map<void *, pair<string, size_t> >::iterator iter = g_AllocMap->begin();
		 iter != g_AllocMap->end(); ++iter)
		{
		const pair<string, size_t> &Pair = iter->second;
		Locvec.push_back(Pair.first);
		Bytesvec.push_back(uint(Pair.second));
		}
	uint N = SIZE(Locvec);
	uint *Order = (uint *) malloc(N*sizeof(uint));
	QuickSortOrderDesc(Bytesvec.data(), N, Order);
	Log("LeakCheck(%s)\n", Msg);
	uint last = UINT_MAX;
	uint64 Sum = 0;
	for (uint k = 0; k < 20u; ++k)
		{
		uint i = Order[k];
		uint bytes = Bytesvec[i];
		Sum += bytes;
		if (bytes > last)
			{
			fprintf(stderr, "ERROR LeakCheck\n");
			exit(1);
			}
		last = bytes;

		if (k < 20)
			{
			const string &Loc = Locvec[i];
			Log("%7u  %s\n", bytes, Loc.c_str());
			}
		}
	free(Order);
	Log("Sum=%s\n", Int64ToStr(Sum));
	}
#endif

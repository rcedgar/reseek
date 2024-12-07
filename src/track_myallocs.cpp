#include "myutils.h"
#include <unordered_map>
#include "sort.h"

#if TRACK_MYALLOCS
extern unordered_map<void *, pair<string, size_t> > *g_AllocMap;

static void LogLo(const char *Msg,
	const vector<string> &Locvec, const vector<uint> &Bytesvec)
	{
	uint N = SIZE(Locvec);
	if (N == 0)
		{
		Log("LogMyAllocs(%s), empty\n", Msg);
		return;
		}
	uint *Order = (uint *) malloc(N*sizeof(uint));
	QuickSortOrderDesc(Bytesvec.data(), N, Order);
	Log("MyAllocs(%s)\n", Msg);
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

void GetMyAllocState(unordered_map<void *, pair<string, size_t> > &AllocMap)
	{
	AllocMap.clear();
	if (g_AllocMap != 0)
		AllocMap = *g_AllocMap;
	}

void CmpMyAllocStates(const char *Msg,
	unordered_map<void *, pair<string, size_t> > &AllocMap1,
	unordered_map<void *, pair<string, size_t> > &AllocMap2)
	{
	vector<string> Locvec1and2;
	vector<uint> Bytesvec1and2;
	vector<string> Locvec2only;
	vector<uint> Bytesvec2only;
	for (unordered_map<void *, pair<string, size_t> >::iterator iter2 = AllocMap2.begin();
		 iter2 != AllocMap2.end(); ++iter2)
		{
		void *Key = iter2->first;
		const pair<string, size_t> &Pair = iter2->second;
		if (AllocMap1.find(Key) == AllocMap1.end())
			{
			Locvec1and2.push_back(Pair.first);
			Bytesvec1and2.push_back(uint(Pair.second));
			}
		else
			{
			Locvec2only.push_back(Pair.first);
			Bytesvec2only.push_back(uint(Pair.second));
			}
		}
	string MsgStr;
	Ps(MsgStr, "%s/1and2", Msg);
	LogLo(MsgStr.c_str(), Locvec1and2, Bytesvec1and2);

	Ps(MsgStr, "%s/2only", Msg);
	LogLo(MsgStr.c_str(), Locvec2only, Bytesvec2only);
	}

void LogMyAllocs(const char *Msg)
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
	LogLo("", Locvec, Bytesvec);
	}
#endif

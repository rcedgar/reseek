#include "myutils.h"
#include "objcounter.h"
#if 0
mutex g_ObjCounter_Lock;
map<string, int> g_ObjCounter_NameToObjCount;
map<string, int> g_ObjCounter_NameToMaxObjCount;

void *ObjCounter_New(size_t n, const string &Name)
	{
	void *p = malloc(n);
	if (p == 0)
		Die("ObjCounter_New");
	g_ObjCounter_Lock.lock();
	map<string, int>::iterator iter = g_ObjCounter_NameToObjCount.find(Name);
	int curr = 0;
	if (iter == g_ObjCounter_NameToObjCount.end())
		{
		g_ObjCounter_NameToObjCount[Name] = 1;
		g_ObjCounter_NameToMaxObjCount[Name] = 1;
		}
	else
		{
		curr = (g_ObjCounter_NameToObjCount[Name] += 1);
		if (curr > g_ObjCounter_NameToMaxObjCount[Name])
			g_ObjCounter_NameToMaxObjCount[Name] = curr;
		}
	g_ObjCounter_Lock.unlock();
	return p;
	}

void ObjCounter_Delete(void *, const string &Name)
	{
	g_ObjCounter_Lock.lock();
	map<string, int>::iterator iter = g_ObjCounter_NameToObjCount.find(Name);
	if (iter == g_ObjCounter_NameToObjCount.end())
		g_ObjCounter_NameToObjCount[Name] = -1;
	else
		g_ObjCounter_NameToObjCount[Name] -= 1;
	g_ObjCounter_Lock.unlock();
	}

void ObjCounter_LogReport()
	{
	g_ObjCounter_Lock.lock();
	Log("\n");
	for (map<string, int>::const_iterator iter = g_ObjCounter_NameToObjCount.begin();
		 iter != g_ObjCounter_NameToObjCount.end(); ++iter)
		{
		const string &Name = iter->first;
		int N = iter->second;
		int Max = g_ObjCounter_NameToMaxObjCount[Name];
		Log("@OBJ %d, %d, %s\n", N, Max, Name.c_str());
		}
	g_ObjCounter_Lock.unlock();
	}
#endif // 0
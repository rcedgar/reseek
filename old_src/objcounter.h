#pragma once

extern mutex g_ObjCounter_Lock;
extern map<string, int> g_ObjCounter_NameToObjCount;

void *ObjCounter_New(size_t n, const string &Name);
void ObjCounter_Delete(void *p, const string &Name);
void ObjCounter_LogReport();

#define OBJCOUNTER(Name)	\
	void *operator new(size_t n) \
		{ \
		void *p = ObjCounter_New(n, #Name); \
		return p; \
		} \
		\
	void operator delete(void *p) \
		{ \
		ObjCounter_Delete(p, #Name); \
		}

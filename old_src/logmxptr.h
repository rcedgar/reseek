#pragma once

template<class t> void LogVecPtr(const string &Name, const char *Fmt,
	const t *VecPtr, uint N)
	{
	Log("\n%s[%u] ", Name.c_str(), N);
	for (uint i = 0; i < N; ++i)
		Log(Fmt, VecPtr[i]);
	Log("\n");
	}

template<class t> void LogMxPtr(const string &Name, const char *Fmt,
	const t *MxPtr, uint N)
	{
	Log("\n%s (%u)\n", Name.c_str(), N);
	for (uint i = 0; i < N; ++i)
		{
		Log("[%3u] ", i);
		for (uint j = 0; j < N; ++j)
			Log(Fmt, MxPtr[N*i + j]);
		Log("\n");
		}
	}

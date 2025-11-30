#include "myutils.h"
#include "dssparams.h"

static vector<string> s_Names;
static vector<string> s_Types;
static vector<void *> s_Ptrs;
static vector<string> s_Defaults;
static vector<bool> s_IsPeakerVarVec;
static map<string, uint> s_NameToIdx;

static bool Init()
	{
	DSSParams::InitParamData();
	return true;
	}
static bool InitDone = Init();

void DSSParams::GetParamsStr(string &ParamsStr)
	{
	ParamsStr.clear();

	const uint FeatureCount = GetFeatureCount();
	asserta(SIZE(m_Features) == FeatureCount);
	asserta(SIZE(m_Weights) == FeatureCount);
	for (uint k = 0; k < FeatureCount; ++k)
		Psa(ParamsStr, "%s=%.3g;",
			FeatureToStr(m_Features[k]),
			m_Weights[k]);

	const uint n = SIZE(s_Names);
	for (uint i = 0; i < n; ++i)
		{
		const string &Name = s_Names[i];
		if (!ParamIsDefault(Name))
			{
			string ValueStr;
			GetParamValueStr(Name, ValueStr);
			ParamsStr += Name + "=" + ValueStr + ";";
			}
		}
	}

bool DSSParams::GetDoMuFilter()
	{
	return GetNeedMuLetters();
	}

bool DSSParams::GetNeedMuLetters()
	{
	if (m_ParaBits == 8)
		return m_Omega8 > 0;
	else if (m_ParaBits == 16)
		return m_Omega16 > 0;
	else if (m_ParaBits == 2)
		{
		asserta(m_Omega8 > 0 && m_Omega16 > 0);
		return true;
		}
	Die("DSSParams::NeedMuLetters");
	return false;
	}

void DSSParams::SetParam(const string &Name, const string &StrValue)
	{
	if (Name == "gap2")
		{
		float Value = StrToFloatf(StrValue);
		m_GapOpen = -Value;
		m_GapExt = -Value/10;
		return;
		}
	if (Name == "open")
		{
		float Value = StrToFloatf(StrValue);
		m_GapOpen = -Value;
		return;
		}
	if (Name == "ext")
		{
		float Value = StrToFloatf(StrValue);
		m_GapExt = -Value;
		return;
		}
	if (Name == "logladd")
		{
		float Value = StrToFloatf(StrValue);
		if (Value < 1)
			Value = 1;
		if (Value > 4)
			Value = 4;
		m_ladd = (float) pow(10, Value);
		return;
		}
	map<string, uint>::const_iterator iter = s_NameToIdx.find(Name);
	if (iter == s_NameToIdx.end())
		Die("DSSParams::SetParam(%s)", Name.c_str());
	uint Idx = iter->second;
	asserta(Idx < SIZE(s_Types));
	const string &Type = s_Types[Idx];
	if (Type == "int")
		{
		int Value = StrToInt(StrValue);
		int *Ptr = (int *) s_Ptrs[Idx];
		*Ptr = Value;
		}
	else if (Type == "uint")
		{
		uint Value = StrToUint(StrValue);
		uint *Ptr = (uint *) s_Ptrs[Idx];
		*Ptr = Value;
		}
	else if (Type == "float")
		{
		float Value = StrToFloatf(StrValue);
		if (Name == "gap2")
			{
			m_GapOpen = -Value;
			m_GapExt = -Value/10;
			}
		else
			{
			float *Ptr = (float *) s_Ptrs[Idx];
			*Ptr = Value;
			}
		}
	else
		Die("DSSParams::SetParam(%s) Type=%s", Name.c_str(), Type.c_str());
	}

int DSSParams::GetIntParam(const string &Name)
	{
	map<string, uint>::const_iterator iter = s_NameToIdx.find(Name);
	asserta(iter != s_NameToIdx.end());
	uint Idx = iter->second;
	void *Ptr = s_Ptrs[Idx];
	int Value = *((int *) Ptr);
	return Value;
	}

uint DSSParams::GetUintParam(const string &Name)
	{
	map<string, uint>::const_iterator iter = s_NameToIdx.find(Name);
	asserta(iter != s_NameToIdx.end());
	uint Idx = iter->second;
	void *Ptr = s_Ptrs[Idx];
	uint Value = *((uint *) Ptr);
	return Value;
	}

float DSSParams::GetFloatParam(const string &Name)
	{
	map<string, uint>::const_iterator iter = s_NameToIdx.find(Name);
	asserta(iter != s_NameToIdx.end());
	uint Idx = iter->second;
	void *Ptr = s_Ptrs[Idx];
	float Value = *((float *) Ptr);
	return Value;
	}

bool DSSParams::ParamIsDefault(const string &Name)
	{
	string ValueStr, DefaultStr;
	GetParamValueStr(Name, ValueStr);
	GetParamDefaultStr(Name, DefaultStr);
	return ValueStr == DefaultStr;
	}

const char *DSSParams::GetParamDefaultStr(const string &Name, string &DefaultStr)
	{
	map<string, uint>::const_iterator iter = s_NameToIdx.find(Name);
	asserta(iter != s_NameToIdx.end());
	const uint Idx = iter->second;
	asserta(Idx < SIZE(s_Defaults));
	DefaultStr = s_Defaults[Idx];
	return DefaultStr.c_str();
	}

const char *DSSParams::GetParamValueStr(const string &Name, string &ValueStr)
	{
	map<string, uint>::const_iterator iter = s_NameToIdx.find(Name);
	asserta(iter != s_NameToIdx.end());
	const uint Idx = iter->second;
	asserta(Idx < SIZE(s_Types));
	const string &Type = s_Types[Idx];
	if (Type == "int")
		{
		int Value = GetIntParam(Name);
		Ps(ValueStr, "%d", Value);
		}
	else if (Type == "uint")
		{
		uint Value = GetUintParam(Name);
		Ps(ValueStr, "%u", Value);
		}
	else if (Type == "float")
		{
		float Value = GetFloatParam(Name);
		Ps(ValueStr, "%.6g", Value);
		}
	else
		Die("Type=%s", Type.c_str());
	return ValueStr.c_str();
	}

void DSSParams::LogParamData()
	{
	const uint n = SIZE(s_Names);
	for (uint i = 0; i < n; ++i)
		{
		const string &Name = s_Names[i].c_str();
		map<string, uint>::const_iterator iter = s_NameToIdx.find(Name);
		asserta(iter != s_NameToIdx.end());
		asserta(iter->second == i);

		const string &Type = s_Types[i];
		string Default = s_Defaults[i];
		if (Type == "int" && Default == "2147483647")
			Default = "INT_MAX";
		if (Type == "uint" && Default == "4294967295")
			Default = "UINT_MAX";
		if (Type == "float" && Default == "3.40282e+38")
			Default = "FLT_MAX";

		Log("%20.20s", Name.c_str());
		Log("  %5.5s", Type.c_str());
		Log("  %10.10s", Default.c_str());

		if (Type == "int")
			{
			int Value = GetIntParam(Name);
			if (Value == StrToInt(s_Defaults[i]))
				Log("  %10.10s", ".");
			else
				Log("  %10d", Value);
			}
		else if (Type == "uint")
			{
			uint Value = GetUintParam(Name);
			if (Value == StrToUint(s_Defaults[i]))
				Log("  %10.10s", ".");
			else
				Log("  %10u", Value);
			}
		else if (Type == "float")
			{
			float Value = GetFloatParam(Name);
			if (feq(Value, StrToFloat(s_Defaults[i])))
				Log("  %10.10s", ".");
			else
				Log("  %10.6g", Value);
			}
		else
			Die("Type=%s", Type.c_str());

		Log("  %c", yon(s_IsPeakerVarVec[i]));
		Log("\n", yon(s_IsPeakerVarVec[i]));
		}
	}

void DSSParams::InitParamData()
	{
#define P(Name)							\
	{									\
	s_Types.push_back("int");			\
	uint Idx = SIZE(s_Names);			\
	s_NameToIdx[#Name] = Idx;			\
	s_Names.push_back(#Name);			\
	s_Ptrs.push_back(&m_##Name);		\
	string Tmp;							\
	Ps(Tmp, "%d", m_##Name);			\
	s_Defaults.push_back(Tmp);			\
	s_IsPeakerVarVec.push_back(false);	\
	}
#include "dssintparams.h"

#define P(Name)							\
	{									\
	s_Types.push_back("uint");			\
	uint Idx = SIZE(s_Names);			\
	s_NameToIdx[#Name] = Idx;			\
	s_Names.push_back(#Name);			\
	s_Ptrs.push_back(&m_##Name);		\
	string Tmp;							\
	Ps(Tmp, "%u", m_##Name);			\
	s_Defaults.push_back(Tmp);			\
	s_IsPeakerVarVec.push_back(false);	\
	}
#include "dssuintparams.h"

#define P(Name)							\
	{									\
	s_Types.push_back("float");			\
	uint Idx = SIZE(s_Names);			\
	s_NameToIdx[#Name] = Idx;			\
	s_Names.push_back(#Name);			\
	s_Ptrs.push_back(&m_##Name);		\
	string Tmp;							\
	Ps(Tmp, "%.6g", m_##Name);			\
	s_Defaults.push_back(Tmp);			\
	s_IsPeakerVarVec.push_back(false);	\
	}
#include "dssfloatparams.h"
	}

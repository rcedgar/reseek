#include "myutils.h"
#include "peaker.h"
#include "sort.h"
#include <time.h>

/***
$res/repeak/misc/spec.txt
=========================
cmd=C:\cygwin64\bin\bash.exe -e z:/a/res/repeak/misc/eval.bash
mindy=0.01
maxdy=1
minh=0.01
sigfig=3
latin=no
var=open        min=0   max=1   delta=0.05      bins=4  init=2
var=ext min=0   max=1   delta=0.05      bins=4  init=0.25
var=topen       min=0   max=1   delta=0.05      bins=4  init=0
var=text        min=0   max=1   delta=0.05      bins=4  init=0

$res/repeak/misc/spec2.txt
==========================
dir=queue/
dy=3200
mindy=50
cool=0.8
sigfig=5
latin=200
var=dpw min=0   max=1   mind=0.01       maxd=0.1        bins=16
var=adsw        min=0   max=1   mind=0.01       maxd=0.2        bins=16
var=lddtw       min=0   max=1   mind=0.01       maxd=0.2        bins=16
var=lddtMw      min=0   max=1   mind=0.01       maxd=0.2        bins=16
***/

FILE *Peaker::m_fTsv = 0;

void Peaker::GetGlobalStr(const string &Name, string &s, const string &Default) const
	{
	return SpecGetStr(m_GlobalSpec, Name, s, Default);
	}

double Peaker::GetGlobalFloat(const string &Name, double Default) const
	{
	return SpecGetFloat(m_GlobalSpec, Name, Default);
	}

uint Peaker::GetGlobalInt(const string &Name, uint Default) const
	{
	return SpecGetInt(m_GlobalSpec, Name, Default);
	}

bool Peaker::GetGlobalBool(const string &Name, bool Default) const
	{
	return SpecGetBool(m_GlobalSpec, Name, Default);
	}

void Peaker::VarSpecGetStr(uint VarIdx, const string &Name,
	string &s, const string &Default) const
	{
	SpecGetStr(GetVarSpec(VarIdx), Name, s, Default);
	}

double Peaker::VarSpecGetFloat(uint VarIdx, const string &Name, double Default) const
	{
	return SpecGetFloat(GetVarSpec(VarIdx), Name, Default);
	}

uint Peaker::VarSpecGetInt(uint VarIdx, const string &Name, uint Default) const
	{
	return SpecGetInt(GetVarSpec(VarIdx), Name, Default);
	}

bool Peaker::VarSpecGetBool(uint VarIdx, const string &Name, bool Default) const
	{
	return SpecGetBool(GetVarSpec(VarIdx), Name, Default);
	}

double Peaker::SpecGetFloat(const string &Spec, const string &Name, double Default)
	{
	string s;
	SpecGetStr(Spec, Name, s, "");
	if (s == "")
		return Default;
	return StrToFloat(s);
	}

void Peaker::SpecGetStr(const string &Spec, const string &Name,
		string &Str, const string &Default)
	{
	vector<string> Fields;
	Split(Spec, Fields, ';');
	const string NameEq = Name + "=";
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		const string &Field = Fields[i];
		if (StartsWith(Field, NameEq))
			{
			vector<string> Fields2;
			Split(Field, Fields2, '=');
			asserta(SIZE(Fields2) == 2);
			Str = Fields2[1];
			return;
			}
		}
	Str = Default;
	}

bool Peaker::SpecGetBool(const string &Spec, const string &Name, bool Default)
	{
	string s;
	SpecGetStr(Spec, Name, s, "");
	if (s == "")
		return Default;
	if (s == "yes")
		return true;
	else if (s == "no")
		return false;
	Die("Peaker::SpecGetBool(%s)", s.c_str());
	return false;
	}

uint Peaker::SpecGetInt(const string &Spec, const string &Name, uint Default)
	{
	string s;
	SpecGetStr(Spec, Name, s, "");
	if (s == "")
		return Default;
	return StrToUint(s);
	}

void Peaker::str2xv(const string &xstr, vector<string> &xv) const
	{
	xv.clear();
	vector<string> Fields;
	Split(xstr, Fields, ';');
	const uint n = SIZE(Fields);
	const uint VarCount = GetVarCount();
	xv.resize(VarCount);
	if (n != VarCount)
		Die("%u/%u str2xv(%s)", n, VarCount, xstr.c_str());
	for (uint k = 0; k < VarCount; ++k)
		{
		const string &NameEqValue = Fields[k];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			Die("Bad name=value %s str2xv(%s)",
				NameEqValue.c_str(), xstr.c_str());
		uint VarIdx = GetVarIdx(Fields2[0]);
		asserta(VarIdx < SIZE(xv));
		asserta(xv[VarIdx] == "");
		xv[VarIdx] = Fields2[1];
		}
	}

void Peaker::xv2str(const vector<string> &xv, string &s) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv) == VarCount);
	s.clear();
	for (uint i = 0; i < VarCount; ++i)
		s += m_VarNames[i] + "=" + xv[i] + ";";
	}

double Peaker::rr(double lo, double hi) const
	{
	const uint M = 18409199;
	double r = (randu32()%(M+1))/double(M);
	double x = lo + r*(hi - lo);
	assert(x >= lo && x <= hi);
	return x;
	}

const char *Peaker::GetVarName(uint VarIdx) const
	{
	asserta(VarIdx < SIZE(m_VarNames));
	return m_VarNames[VarIdx].c_str();
	}

double Peaker::GetBinWidth(uint VarIdx) const
	{
	if (VarIsConstant(VarIdx))
		return 0;
	uint BinCount = GetGlobalInt("latin", 0);
	if (BinCount == 0)
		return 0;
	double Min = VarSpecGetFloat(VarIdx, "min", true);
	double Max = VarSpecGetFloat(VarIdx, "max", true);
	if (Min == DBL_MAX || Max == DBL_MAX)
		return 0;

	asserta(Max > Min);
	double BinWidth = (Max - Min)/BinCount;
	return BinWidth;
	}

double Peaker::GetLatinValueByBinIdx(uint VarIdx, uint BinIdx, uint BinCount) const
	{
	asserta(BinCount > 0);
	double Min = VarSpecGetFloat(VarIdx, "min", DBL_MAX);
	double Max = VarSpecGetFloat(VarIdx, "max", DBL_MAX);
	asserta(Min != DBL_MAX && Max != DBL_MAX);
	double BinWidth = (Max - Min)/BinCount;
	asserta(feq(BinWidth*BinCount, Max - Min));
	double BinLo = Min + BinWidth*BinIdx;
	double rf = randf(1);
	double Value = BinLo + rf*BinWidth;
	asserta(Value >= Min && Value <= Max);
	return Value;
	}

void Peaker::GetLatinHypercube(vector<vector<string > > &xvs)
	{
	xvs.clear();
	uint BinCount = GetGlobalInt("latin", 0);
	if (BinCount == 0)
		return;

	const uint VarCount = GetVarCount();
	vector<vector<uint> > IdxMx;
	GetLatinHypercubeIdxs(VarCount, BinCount, IdxMx);
	asserta(SIZE(IdxMx) == BinCount);

	xvs.resize(BinCount);
	for (uint k = 0; k < BinCount; ++k)
		{
		vector<string> &xv = xvs[k];
		for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
			{
			uint BinIdx = IdxMx[k][VarIdx];
			double Value = GetLatinValueByBinIdx(VarIdx, BinIdx, BinCount);
			string ValueStr;
			VarFloatToStr(VarIdx, Value, ValueStr);
			xv.push_back(ValueStr);
			}
		}
	}

const string &Peaker::GetVarSpec(uint VarIdx) const
	{
	asserta(VarIdx < SIZE(m_VarSpecs));
	return m_VarSpecs[VarIdx];
	}

bool Peaker::VarIsConstant(uint VarIdx) const
	{
	string yes;
	VarSpecGetStr(VarIdx, "constant", yes, "no");
	return yes == "yes";
	}

uint Peaker::GetVarIdx(const string &Name, bool FailOk) const
	{
	for (uint Idx = 0; Idx < SIZE(m_VarNames); ++Idx)
		{
		if (Name == m_VarNames[Idx])
			return Idx;
		}
	if (FailOk)
		return UINT_MAX;
	Die("GetVarIdx(%s)", Name.c_str());
	return 0;
	}

void Peaker::InitRates()
	{
	m_Rates.clear();
	const uint VarCount = GetVarCount();
	m_Rates.resize(VarCount, MED_RATE);
	}

double Peaker::Calc(const vector<string> &xv)
	{
	asserta(m_EvalFunc != 0);
	double y = (*m_EvalFunc)(xv);
	return y;
	}

void Peaker::GetPeakerPathStr(string &s) const
	{
	s.clear();
	string ParentStr;
	if (m_Parent != 0)
		{
		m_Parent->GetPeakerPathStr(ParentStr);
		if (ParentStr == "/")
			ParentStr = "";
		}

	string WhoAmI;
	s = ParentStr + "/" + WhoAmI;
	}

void Peaker::AppendChildResult(const vector<string> &xv, double y, 
	const string &childname, const string &why)
	{
	m_xvs.push_back(xv);
	m_ys.push_back(y);

	string child_why = childname + ":" + why;
	m_whys.push_back(child_why);

	string improved;
	if (y != DBL_MAX)
		{
		if (m_Best_y == DBL_MAX || y > m_Best_y)
			{
			improved = "child>>";
			m_Best_y = y;
			m_Best_xv = xv;
			m_LastImprovedTime = time(0);
			}
		}
	const uint n = SIZE(m_xvs);
	asserta(SIZE(m_whys) == n);
	asserta(SIZE(m_ys) == n);
	if (m_fTsv != 0)
		{
		string tmp;
		xv2str(xv, tmp);
		fprintf(m_fTsv, "%.6g\t%s%s\t%s\n",
			y, improved.c_str(), why.c_str(), tmp.c_str());
		fflush(m_fTsv);
		}
	}

uint Peaker::Find_xv(const vector<string> &xv) const
	{
	const uint N = SIZE(m_xvs);
	asserta(SIZE(m_ys) == N);
	for (uint Idx = 0; Idx < N; ++Idx)
		{
		if (m_xvs[Idx] == xv)
			return Idx;
		}
	return UINT_MAX;
	}

void Peaker::NormalizeWeights(const vector<string> &xv,
	vector<string> &Normalized_xv) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv) == VarCount);
	vector<double> Weights;
	double SumWeight = 0;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		if (!VarSpecGetBool(VarIdx, "weight", false))
			continue;
		double Weight = StrToFloat(xv[VarIdx]);
		Weights.push_back(Weight);
		SumWeight += Weight;
		}
	Normalized_xv = xv;
	if (SIZE(Weights) == 0)
		return;
	asserta(SumWeight > 0);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		if (!VarSpecGetBool(VarIdx, "weight", false))
			continue;
		uint SigFig = VarSpecGetInt(VarIdx, "sigfig", 2);
		double Weight = StrToFloat(xv[VarIdx]);
		Weight /= SumWeight;
		string TmpStr;
		GetRoundedStr(Weight, SigFig, TmpStr);
		string NewStr;
		NormalizeVarStr(VarIdx, TmpStr, NewStr);
		Normalized_xv[VarIdx] = NewStr;
		}
	}

double Peaker::Evaluate(const vector<string> &axv, const string &why)
	{
	vector<string> xv;
	NormalizeWeights(axv, xv);
	const double Saved_Best_y = m_Best_y;
	uint Idx = Find_xv(xv);
	double y = DBL_MAX;
	if (Idx != UINT_MAX)
		{
		++m_EvaluateCacheHits;
		asserta(Idx < SIZE(m_ys));
		y = m_ys[Idx];
		}
	else
		y = Calc(xv);

	m_xvs.push_back(xv);
	m_ys.push_back(y);
	m_whys.push_back(why);

	string Path;
	GetPeakerPathStr(Path);

	if (y != DBL_MAX)
		{
		if (m_Best_y == DBL_MAX || y > m_Best_y)
			{
			m_Best_y = y;
			m_Best_xv = xv;
			m_LastImprovedTime = time(0);
			m_LastImprovedEvalCount = SIZE(m_ys);
			ProgressLog(">>>%.2g [%.5g] %s %s\n",
				y - Saved_Best_y, y, Path.c_str(), why.c_str());
			}
		}

	double dy = 0;
	if (Saved_Best_y != DBL_MAX)
		dy = y - Saved_Best_y;
	string tmp;
	xv2str(xv, tmp);
	ProgressLog("[%.6g](%+.3g) %s %s %s\n",
		m_Best_y, dy, Path.c_str(), why.c_str(), tmp.c_str());
	string improved;
	if (dy > 0)
		improved = ">>";
	if (m_fTsv != 0 && y != DBL_MAX)
		{
		string tmp;
		xv2str(xv, tmp);
		fprintf(m_fTsv, "%.6g\t%s%s:%s\t%s\n",
			y, improved.c_str(), Path.c_str(), why.c_str(), tmp.c_str());
		fflush(m_fTsv);
		}
	return y;
	}

void Peaker::Init(const vector<string> &SpecLines, PTR_EVAL_FUNC EF)
	{
	m_GlobalSpec.clear();
	m_VarSpecs.clear();
	m_EvalFunc = EF;
	vector<string> Fields;
	const uint LineCount = SIZE(SpecLines);
	for (uint LineIdx = 0; LineIdx < LineCount; ++LineIdx)
		{
		const string &Line = SpecLines[LineIdx];
		if (Line.empty() || Line[0] == '#')
			continue;
		if (StartsWith(Line, "var="))
			{
			string Name;
			SpecGetStr(Line, "var", Name, "");
			asserta(!Name.empty());
			m_VarNames.push_back(Name);
			m_VarSpecs.push_back(Line);
			}
		else
			{
			if (StartsWith(Line, "init="))
				{
				asserta(m_InitParams.empty());
				m_InitParams = Line.substr(5);
				}
			else
				{
				asserta(EndsWith(Line, ";"));
				m_GlobalSpec += Line;
				}
			}
		}
	}

void Peaker::LogState() const
	{
	Die("TODO");
	}

double Peaker::GetRoundedValue(double x, uint SigFig)
	{
	string Fmt;
	Ps(Fmt, "%%.%uE", SigFig-1);
	string Tmp;
	Ps(Tmp, Fmt.c_str(), x);
	double Rounded_x = atof(Tmp.c_str());
	return Rounded_x;
	}

void Peaker::GetRoundedStr(double x, uint SigFig, string &Str)
	{
	asserta(SigFig > 1);
	string Fmt;
	Ps(Fmt, "%%.%uE", SigFig-1);
	Ps(Str, Fmt.c_str(), x);
	}

void Peaker::RunLatin()
	{
	uint BinCount = GetGlobalInt("latin", 0);
	if (BinCount == 0)
		{
		ProgressLog("RunLatin(), no bins\n");
		return;
		}

	const uint VarCount = GetVarCount();
	vector<vector<string> > xvs;
	GetLatinHypercube(xvs);
	const uint n = SIZE(xvs);
	asserta(n > 0);
#if	1
	{
	const uint m = GetVarCount();
	Log("Latin=%u\n", n);
	for (uint i = 0; i < n; ++i)
		{
		const vector<string> &xv = xvs[i];
		string tmp;
		xv2str(xv, tmp);
		Log("[%5u] %s", i, tmp.c_str());
		Log("\n");
		}
	}
#endif
	for (uint i = 0; i < n; ++i)
		{
		string why;
		Ps(why, "latin%u/%u", i+1, n);
		Evaluate(xvs[i], why);
		}
	}

void Peaker::GetLatinHypercubeIdxs(uint VarCount, uint BinCount,
	vector<vector<uint> > &IdxMx)
	{
	IdxMx.clear();
	IdxMx.resize(BinCount);

	for (uint BinIdx = 0; BinIdx < BinCount; ++BinIdx)
		IdxMx[BinIdx].resize(VarCount, BinIdx);

	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		vector<uint> BinIdxs;
		for (uint BinIdx = 0; BinIdx < BinCount; ++BinIdx)
			BinIdxs.push_back(BinIdx);
		Shuffle(BinIdxs);
		for (uint k = 0; k < BinCount; ++k)
			IdxMx[k][VarIdx] = BinIdxs[k];
		}
	}

void Peaker::WriteStatusPage(FILE *f) const
	{
	if (f == 0)
		return;
	}

void Peaker::LogSpec() const
	{
	const uint n = SIZE(m_SpecLines);
	Log("LogSpec %u lines %s\n", n, m_SpecFN.c_str());
	for (uint i = 0; i < n; ++i)
		Log("%s\n", m_SpecLines[i].c_str());
	}

void Peaker::RunNestedLatin(uint TopN)
	{
	ProgressLog("RunNestedLatin(%u)\n", TopN);
	RunLatin();

	const vector<double> &ys = m_ys;
	const vector<vector<string> > &xvs = m_xvs;
	const uint NY = SIZE(ys);
	const uint N = min(TopN, NY);
	if (N == 0)
		return;
	const uint VarCount = GetVarCount();
	vector<uint> Order(NY);
	QuickSortOrderDesc(ys.data(), NY, Order.data());
	Log("NESTED_LATIN_TOPS ");
	for (uint k = 0; k < N; ++k)
		{
		uint i = Order[k];
		Log(" y=%.3g", ys[i]);
		}
	Log("\n");
	const double yseed0 = ys[Order[0]];

	for (uint k = 0; k < N; ++k)
		{
		uint i = Order[k];
		asserta(i < SIZE(m_xvs));
		const vector<string> &xv = m_xvs[i];
		double yseed = ys[i];
		asserta(SIZE(xv) == VarCount);
		Log("RUN_NESTED_LATIN %u/%u y=%.3g\n", k+1, N, ys[i]);

		string nest;
		Ps(nest, "NestedLatin%u", k+1);
		Peaker P(this, nest);
		P.m_GlobalSpec = m_GlobalSpec;
		P.m_EvalFunc = m_EvalFunc;
		P.m_fTsv = m_fTsv;
		for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
			{
			const string &VarStr = xv[VarIdx];
			double BinWidth = GetBinWidth(VarIdx);
			string NewSpec = GetVarSpec(VarIdx);
			double CenterValue = VarStrToFloat(VarIdx, VarStr);
			double NewMin = CenterValue - BinWidth;
			if (NewMin < 0)
				NewMin = 0;
			double NewMax = CenterValue + BinWidth;
			string NewMinStr, NewMaxStr;
			VarFloatToStr(VarIdx, NewMin, NewMinStr);
			VarFloatToStr(VarIdx, NewMax, NewMaxStr);
			asserta(NewMinStr != NewMaxStr);
			SpecReplaceStr(NewSpec, "min", NewMinStr);
			SpecReplaceStr(NewSpec, "max", NewMaxStr);
			P.m_VarSpecs.push_back(NewSpec);
			}

		Log("\n");
		Log("NESTED_LATIN\n");
		P.LogSpec();
		P.RunLatin();
		const vector<vector<string> > &Child_xvs = P.m_xvs;
		const vector<double> &Child_ys = P.m_ys;
		const uint J = SIZE(Child_xvs);
		asserta(SIZE(P.m_whys) == J);
		asserta(SIZE(Child_ys) == J);
		Log("NESTED_LATIN_DONE %u/%u yseed0 %.3g yseed=%.3g best=%.3g J=%u\n",
			k+1, N, yseed0, yseed, P.m_Best_y, J);
		for (uint j = 0; j < J; ++j)
			{
			const vector<string> &Child_xv = Child_xvs[j];
			double y = Child_ys[j];
			const string &Child_why = P.m_whys[j];
			AppendChildResult(Child_xv, y, nest, Child_why);
			}
		}
	}

double Peaker::VarStrToFloat(uint VarIdx, const string &ValueStr) const
	{
	return StrToFloat(ValueStr);
	}

void Peaker::VarFloatToStr(uint VarIdx, double Value, string &Str) const
	{
	asserta(Value != DBL_MAX);
	asserta(!isnan(Value));
	uint SigFig = VarSpecGetInt(VarIdx, "sigfig", 2);

	string Fmt;
	Ps(Fmt, "%%.%ug", SigFig);

	string TmpStr;
	Ps(TmpStr, Fmt.c_str(), Value);

	NormalizeVarStr(VarIdx, TmpStr, Str);
	}

void Peaker::SpecReplaceStr(string &Spec, const string &Name,
	const string &NewValue)
	{
	vector<string> Fields;
	Split(Spec, Fields, ';');
	Spec.clear();
	const string NameEq = Name + "=";
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		const string &Field = Fields[i];
		if (StartsWith(Field, NameEq))
			{
			vector<string> Fields2;
			Split(Field, Fields2, '=');
			asserta(SIZE(Fields2) == 2);
			const string &Name = Fields2[0];
			Spec += NameEq + NewValue + ";";
			}
		else
			Spec += Fields[i];
		}
	}

uint Peaker::GetSigFig(const string &EStr)
	{
	Die("TODO");
	return 0;
	}

void Peaker::ParseEStr(const string &EStr, string &Mantissa, string &Exponent)
	{
	// 1E1
	asserta(SIZE(EStr) >= 3);
	Mantissa.clear();
	Exponent.clear();
	const uint n = SIZE(EStr);
	uint i = 0;
	while (i < n)
		{
		char c = EStr[i++];
		if (i == 2)
			{
			asserta(c == '.');
			continue;
			}
		if (c == 'E' || c == 'e')
			{
			asserta(!Mantissa.empty());
			break;
			}
		asserta(isdigit(c));
		Mantissa += c;
		}
	if (Mantissa[0] == '0')
		{
		for (uint i = 1; i < SIZE(Mantissa); ++i)
			asserta(Mantissa[i] == '0');
		}
	while (i < n)
		{
		char c = EStr[i++];
		asserta(c == '-' || c == '+' || isdigit(c));
		Exponent += c;
		}
	asserta(!Exponent.empty());
	}

// All chars are '0'
bool Peaker::AllZeros(const string &s)
	{
	const uint n = SIZE(s);
	if (n == 0)
		return false;
	for (uint i = 0; i < n; ++i)
		if (s[i] != '0')
			return false;
	return true;
	}

// All chars are '9'
bool Peaker::AllNines(const string &s)
	{
	const uint n = SIZE(s);
	if (n == 0)
		return false;
	for (uint i = 0; i < n; ++i)
		if (s[i] != '9')
			return false;
	return true;
	}

// First char is '1' rest (if any) are '0'
bool Peaker::OneZeros(const string &s)
	{
	const uint n = SIZE(s);
	if (n == 0)
		return false;
	if (s[0] != '1')
		return false;
	for (uint i = 1; i < n; ++i)
		if (s[i] != '0')
			return false;
	return true;
	}

void Peaker::IncFloat(const string &OldStr, bool Plus, string &NewStr)
	{
	NewStr.clear();
	double OldValue = StrToFloat(OldStr);

	string mantissa_str;
	string exponent_str;
	ParseEStr(OldStr, mantissa_str, exponent_str);

	uint imantissa = StrToUint(mantissa_str);
	const size_t mantissa_str_size = mantissa_str.size();
	uint new_imantissa = (Plus ? imantissa + 1 : imantissa - 1);
	string new_imantissa_str;
	Ps(new_imantissa_str, "%u", new_imantissa);
	size_t new_imantissa_str_size = SIZE(new_imantissa_str);
	if (new_imantissa_str_size == SIZE(mantissa_str))
		{
		string new_mantissa_str;
		for (uint i = 0; i < new_imantissa_str_size; ++i)
			{
			new_mantissa_str += new_imantissa_str[i];
			if (i == 0)
				new_mantissa_str += '.';
			}
		NewStr = new_mantissa_str + "E" + exponent_str;
		}
	else if (Plus && AllNines(mantissa_str))
		{
		int new_iexponent = atoi(exponent_str.c_str()) + 1;
		string new_exponent_str;
		Ps(new_exponent_str, "%d", new_iexponent);
		string new_mantissa_str;
		for (uint i = 0; i < mantissa_str_size; ++i)
			{
			if (i == 0)
				new_mantissa_str += "1.";
			else
				new_mantissa_str += '0';
			}
		NewStr = new_mantissa_str + "E" + new_exponent_str;
		}
	else if (!Plus && OneZeros(mantissa_str))
		{
		int new_iexponent = atoi(exponent_str.c_str()) - 1;
		string new_exponent_str;
		Ps(new_exponent_str, "%d", new_iexponent);
		string new_mantissa_str;
		for (uint i = 0; i < mantissa_str_size; ++i)
			{
			if (i == 0)
				new_mantissa_str += "9.";
			else
				new_mantissa_str += '9';
			NewStr = new_mantissa_str + "E" + new_exponent_str;
			}
		}
	// Cannot inc 0.0E0
	else if (AllZeros(mantissa_str))
		NewStr = OldStr;
	else
		asserta(false);
	asserta(!NewStr.empty());
	if (Plus)
		asserta(StrToFloat(NewStr) > StrToFloat(OldStr));
	else
		asserta(StrToFloat(NewStr) < StrToFloat(OldStr));
	}

static void Test1(const string &s, bool Plus)
	{
	string news;
	Peaker::IncFloat(s, Plus, news);
	ProgressLog("%8.8s  %c  %s\n", s.c_str(), pom(Plus), news.c_str());
	}

void cmd_test()
	{
	Test1("9.9E1", true);
	Test1("9.99E2", true);
	Test1("9.999E2", true);
	Test1("9.999E40", true);
	Test1("1.0E1", false);
	Test1("1.00E2", false);
	Test1("1.000E3", false);
	Test1("1.0000E40", false);

	for (uint i = 0; i < 100; ++i)
		{
		uint SigFig = randu32()%3 + 2;
		double mant = randf(1);
		double e = randf(20);
		bool sign = (randu32()%2 == 0);
		double p = (sign ? pow(10, e) : pow(10, -e));
		double x = mant*p;
		string OldStr;
		Peaker::GetRoundedStr(x, SigFig, OldStr);
		bool Plus = (randu32()%2 == 0);
		string NewStr;
		Peaker::IncFloat(OldStr, Plus, NewStr);
		ProgressLog("%10.10s  %c  %s\n",
			OldStr.c_str(), pom(Plus), NewStr.c_str());
		}
	}

#include "myutils.h"
#include "peaker.h"
#include "sort.h"
#include <time.h>

FILE *Peaker::m_fTsv = 0;

void Peaker::GetGlobalSpec(const vector<string> &SpecLines, string &GlobalSpec)
	{
	GlobalSpec.clear();
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		{
		const string &Line = SpecLines[i];
		Log("%s\n", Line.c_str());
		if (Line.empty())
			continue;
		if (StartsWith(Line, "#"))
			continue;
		if (StartsWith(Line, "var="))
			continue;
		GlobalSpec += Line;
		}
	}

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
			if (SIZE(Fields2) != 2)
				Die("expected name=value '%s'", Field.c_str());
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

// err_ok=true for resuming, tsv may be prematurely truncated
bool Peaker::xss2xv(const string &xstr, vector<string> &xv, bool err_ok) const
	{
	xv.clear();
	vector<string> Fields;
	Split(xstr, Fields, ';');
	const uint n = SIZE(Fields);
	const uint VarCount = GetVarCount();
	xv.resize(VarCount);
	if (n != VarCount)
		{
		if (err_ok) return false;
		Die("%u/%u xss2xv(%s)", n, VarCount, xstr.c_str());
		}
	for (uint k = 0; k < VarCount; ++k)
		{
		const string &NameEqValue = Fields[k];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			{
			if (err_ok) return false;
			Die("Bad name=value %s xss2xv(%s)",
				NameEqValue.c_str(), xstr.c_str());
			}
		uint VarIdx = GetVarIdx(Fields2[0]);
		asserta(VarIdx < SIZE(xv));
		asserta(xv[VarIdx] == "");
		xv[VarIdx] = Fields2[1];
		}
	return true;
	}

const char *Peaker::xv2xss(const vector<string> &xv, string &s) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv) == VarCount);
	s.clear();
	for (uint i = 0; i < VarCount; ++i)
		s += m_VarNames[i] + "=" + xv[i] + ";";
	return s.c_str();
	}

double Peaker::rr(double lo, double hi)
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
	if (VarIsConstant(VarIdx))
		{
		double Value = VarSpecGetFloat(VarIdx, "constant", DBL_MAX);
		asserta(Value != DBL_MAX);
		return Value;
		}

	if (VarIsInt(VarIdx))
		{
		int Min = VarSpecGetInt(VarIdx, "min", INT_MAX);
		int Max = VarSpecGetInt(VarIdx, "max", INT_MAX);
		asserta(Min != INT_MAX);
		asserta(Max != INT_MAX);
		double BinWidth = double(Max - Min)/BinCount;
		double Value = round(Min + BinWidth*BinIdx);
		asserta(Value >= Min && Value <= Max);
		return Value;
		}

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

void Peaker::GetLatinHypercube(uint BinCount,
	vector<vector<string > > &xvs) const
	{
	xvs.clear();
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
	string value;
	VarSpecGetStr(VarIdx, "constant", value, "");
	return value != "";
	}

bool Peaker::VarIsInt(uint VarIdx) const
	{
	string Type;
	VarSpecGetStr(VarIdx, "type", Type, "float");
	return Type == "int";
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
	string s;
	GetGlobalStr("rates", s, "1.1");
	m_GlobalVarRateFactorIdx = 0;
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
	if (m_Parent == 0)	
		{
		s = m_Name;
		return;
		}
	s = m_Parent->m_Name + "." + m_Name;
	}

void Peaker::AppendChildResults(const Peaker &Child)
	{
	const string &childname = Child.m_Name;
	m_ChildNames.push_back(childname);
	const uint ChildEvalCount = SIZE(Child.m_xvs);
	asserta(SIZE(Child.m_ys) == ChildEvalCount);
	asserta(SIZE(Child.m_whys) == ChildEvalCount);

	for (uint Idx = 0; Idx < ChildEvalCount; ++Idx)
		{
		const vector<string> &xv = Child.m_xvs[Idx];
		const double y  = Child.m_ys[Idx];
		const string &childwhy  = Child.m_whys[Idx];
		const string &why = childname + ":" + childwhy;
		AppendResult(xv, y, why);
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

void Peaker::WriteFinalResults(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "\n_____________________________________________\n");
	string best_xss;
	xv2xss(m_Best_xv, best_xss);
	fprintf(f, "FINAL %s [%.6g] %s\n",
		m_Name.c_str(), m_Best_y, best_xss.c_str());
	fprintf(f, "\n_____________________________________________\n");
	WriteFinalPeak(f);
	fflush(f);
	}

void Peaker::AppendResult(const vector<string> &xv, double y,
	const string &why)
	{
	if (y == DBL_MAX)
		return;

	const double Saved_Best_y = m_Best_y;
	m_xvs.push_back(xv);
	m_ys.push_back(y);
	m_whys.push_back(why);

	string Path;
	GetPeakerPathStr(Path);
	string desc = (Path.empty() ? why : Path + ":" + why);

	double dy = 0;
	if (m_Best_y != DBL_MAX)
		dy = y - m_Best_y;
	if (m_Best_y == DBL_MAX || y > m_Best_y)
		{
		m_Best_y = y;
		m_Best_xv = xv;
		m_Best_desc = desc;
		m_LastImprovedTime = time(0);
		m_LastImprovedEvalCount = SIZE(m_ys);

		m_Best_ys.push_back(y);
		m_Best_xvs.push_back(xv);
		m_Best_descs.push_back(desc);
		}

	string xss;
	xv2xss(xv, xss);
	if (dy > 0)
		ProgressPrefixLog("\n");
	ProgressPrefixLog("%s%.2g%%[%.6g] %s /%.2f/ %s\n",
		(dy > 0 ? ">>>" : ""),
		GetPct(dy, m_Best_y),
		m_Best_y,
		desc.c_str(),
		GetGlobalRateFactor(),
		xss.c_str());
	if (dy > 0)
		ProgressPrefixLog("\n");

	if (m_fTsv != 0)
		{
		fprintf(m_fTsv, "%.6g", y);
		fprintf(m_fTsv, "\t%s", why.c_str());
		fprintf(m_fTsv, "\t%.2g", dy);
		fprintf(m_fTsv, "\t%s", desc.c_str());
		fprintf(m_fTsv, "\t%s", xss.c_str());
		fprintf(m_fTsv, "\n");
		fflush(m_fTsv);
		}
	}

double Peaker::Evaluate(const vector<string> &axv, const string &awhy)
	{
	string why = awhy;
	vector<string> xv;
	NormalizeWeights(axv, xv);
	uint Idx = Find_xv(xv);
	double y = DBL_MAX;
	if (Idx != UINT_MAX)
		{
		++m_EvaluateCacheHits;
		asserta(Idx < SIZE(m_ys));
		y = m_ys[Idx];
		why += ":cache";
		}
	else
		y = Calc(xv);
	AppendResult(xv, y, why);
	return y;
	}

double Peaker::Evaluate(const string &varstr, const string &why)
	{
	vector<string> flds;
	Split(varstr, flds, ';');
	vector<string> xv;
	for (size_t i = 0; i < flds.size(); ++i)
		{
		vector<string> flds2;
		Split(flds[i], flds2, '=');
		asserta(flds2.size() == 2);
		xv.push_back(flds2[1]);
		}
	return Evaluate(xv, why);
	}

void Peaker::GetVarNames(
	const vector<string> &SpecLines,
	vector<string> &VarNames)
	{
	Peaker Ptmp(0, "Ptmp");
	Ptmp.Init(SpecLines, 0);
	VarNames = Ptmp.m_VarNames;
	}

void Peaker::Init(const vector<string> &SpecLines, PTR_EVAL_FUNC EF)
	{
	m_GlobalSpec.clear();
	m_VarNames.clear();
	m_VarSpecs.clear();
	m_EvalFunc = EF;
	m_SpecLines = SpecLines;

	const uint n = SIZE(SpecLines);
	for (uint i = 0; i < n; ++i)
		{
		const string &Line = SpecLines[i];
		if (Line.empty() || StartsWith(Line, "#"))
			continue;
//		asserta(EndsWith(Line, ";"));
		if (StartsWith(Line, "var="))
			m_VarSpecs.push_back(Line);
		else
			m_GlobalSpec += Line;
		}

	vector<string> Fields;
	const uint VarCount = SIZE(m_VarSpecs);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const string &VarSpec = m_VarSpecs[VarIdx];
		asserta(StartsWith(VarSpec, "var="));
		string Name;
		SpecGetStr(VarSpec, "var", Name, "");
		asserta(!Name.empty());
		m_VarNames.push_back(Name);
		}
	}

void Peaker::LogState() const
	{
	Die("TODO");
	}

double Peaker::GetRoundedValue(double x, uint SigFig)
	{
	if (SigFig == 0)
		return round(x);
	string Fmt;
	Ps(Fmt, "%%.%uE", SigFig-1);
	string Tmp;
	Ps(Tmp, Fmt.c_str(), x);
	double Rounded_x = atof(Tmp.c_str());
	return Rounded_x;
	}

void Peaker::GetRoundedStr(double x, uint SigFig, string &Str)
	{
	if (SigFig == 0)
		{
		Ps(Str, "%d", int(round(x)));
		return;
		}

	asserta(SigFig > 1);
	string Fmt;
	Ps(Fmt, "%%.%uE", SigFig-1);
	Ps(Str, Fmt.c_str(), x);
	}

bool Peaker::RunInits(bool Required)
	{
	const uint VarCount = GetVarCount();
	vector<string> xv;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		string s;
		VarSpecGetStr(VarIdx, "init", s, "");
		if (s == "")
			return false;
		xv.push_back(s);
		}
	Evaluate(xv, "Inits");
	return true;
	}

void Peaker::RunLatin(uint BinCount)
	{
	const uint VarCount = GetVarCount();
	vector<vector<string> > xvs;
	GetLatinHypercube(BinCount, xvs);
	const uint n = SIZE(xvs);
	asserta(n > 0);
	const uint m = GetVarCount();

	Log("Latin=%u\n", n);
	for (uint i = 0; i < n; ++i)
		{
		const vector<string> &xv = xvs[i];
		string tmp;
		xv2xss(xv, tmp);
		Log("[%5u] %s", i, tmp.c_str());
		Log("\n");
		}

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

Peaker *Peaker::MakeChild(const string &Name) const
	{
	Peaker &Child = *new Peaker(this, Name);
	Child.m_GlobalSpec = m_GlobalSpec;
	Child.m_EvalFunc = m_EvalFunc;
	Child.m_VarNames = m_VarNames;
	Child.m_VarSpecs = m_VarSpecs;
	Child.m_fTsv = m_fTsv;
	return &Child;
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

void Peaker::IncFloat(const string &OldStr, bool Plus, string &NewStr, uint SigFig)
	{
	if (SigFig == 0)
		{
		int OldValue = int(round(StrToFloat(OldStr)));
		int NewValue = (Plus ? OldValue + 1 : OldValue - 1);
		Ps(NewStr, "%d", NewValue);
		return;
		}

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

void Peaker::xv2values(const vector<string> &xv, vector<double> &Values) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv) == VarCount);
	Values.clear();
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		Values.push_back(VarStrToFloat(VarIdx, xv[VarIdx]));
	}

void Peaker::GetTopEvalIdxs_mindy(const uint N, double mindy,
	vector<uint> &Idxs) const
	{
	Idxs.clear();
	const uint Count = SIZE(m_ys);
	uint *Order = myalloc(uint, Count);

	QuickSortOrderDesc(m_ys.data(), Count, Order);
	ProgressPrefixLog("GetTopEvalIdxs(%u)\n", N);
	double Last_y = FLT_MAX;
	for (uint i = 0; i < Count; ++i)
		{
		uint Idx = Order[i];
		if (SIZE(Idxs) >= N)
			break;
		double y = m_ys[Idx];
		if (Last_y != FLT_MAX && Last_y - y < mindy)
			continue;
		Idxs.push_back(Idx);
		Last_y = y;
		ProgressPrefixLog("Top [%2u]  %.6g\n", i+1, y);
		}

	myfree(Order);
	}

void Peaker::GetTopEvalIdxs(const uint N, vector<uint> &Idxs) const
	{
	Idxs.clear();
	const uint Count = SIZE(m_ys);
	uint *Order = myalloc(uint, Count);

	QuickSortOrderDesc(m_ys.data(), Count, Order);
	uint n = min(N, Count);
	ProgressPrefixLog("GetTopEvalIdxs(%u)\n", n);
	for (uint i = 0; i < n; ++i)
		{
		uint Idx = Order[i];
		Idxs.push_back(Idx);
		ProgressPrefixLog("Top [%2u]  %.3g\n", i+1, m_ys[Idx]);
		}

	myfree(Order);
	}

double Peaker::GetEuclideanDist(const vector<string> &xv1,
	const vector<string> &xv2) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv1) == VarCount && SIZE(xv2) == VarCount);
	double Sum2 = 0;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		double x1 = VarStrToFloat(VarIdx, xv1[VarIdx]);
		double x2 = VarStrToFloat(VarIdx, xv2[VarIdx]);
		double d = x1 - x2;
		Sum2 += d*d;
		}
	return sqrt(Sum2);
	}

double Peaker::GetEuclideanDist(const vector<double> &xv1,
	const vector<double> &xv2) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv1) == VarCount && SIZE(xv2) == VarCount);
	double Sum2 = 0;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		double d = xv1[VarIdx] - xv2[VarIdx];
		Sum2 += d*d;
		}
	return sqrt(Sum2);
	}

bool Peaker::GetNearestNeighbor(const vector<string> &xv,
	uint VarIdx, bool Plus, vector<string> &Neighbor_xv, double &y) const
	{
	y = -1;
	Neighbor_xv.clear();
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv) == VarCount);
	asserta(VarIdx < VarCount);
	const uint N = SIZE(m_xvs);
	double MinDist = DBL_MAX;
	double Value = StrToFloat(xv[VarIdx]);
	for (uint i = 0; i < N; ++i)
		{
		const vector<string> &tmp_xv = m_xvs[i];
		if (tmp_xv == xv)
			continue;
		double tmp_Value = StrToFloat(tmp_xv[VarIdx]);
		if (Plus)
			{
			if (tmp_Value <= Value)
				continue;
			}
		else
			{
			if (tmp_Value >= Value)
				continue;
			}
		double Dist = GetEuclideanDist(xv, tmp_xv);
		if (Dist < MinDist)
			{
			Neighbor_xv = tmp_xv;
			MinDist = Dist;
			y = m_ys[i];
			}
		}
	if (MinDist == DBL_MAX)
		return false;
	return true;
	}

void Peaker::WriteFinalPeak(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "Final peak\n");
	const uint VarCount = GetVarCount();
	vector<string> Neighbor;
	string xss;
	double y;
	vector<string> MinValueStrs(VarCount);
	vector<string> MaxValueStrs(VarCount);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		fprintf(f, "\n");
		fprintf(f, "%s\n", GetVarName(VarIdx));
		for (int iPlus = 0; iPlus <= 1; ++iPlus)
			{
			bool Plus = (iPlus == 0);
			bool Ok = GetNearestNeighbor(m_Best_xv, VarIdx, Plus, Neighbor, y);
			fprintf(f, "%c", pom(Plus));
			if (!Ok)
				{
				fprintf(f, "  (none)\n");
				continue;
				}
			for (uint VarIdx2 = 0; VarIdx2 < VarCount; ++VarIdx2)
				{
				double Value = StrToFloat(Neighbor[VarIdx2]);
				if (MinValueStrs[VarIdx2].empty() || Value < StrToFloat(MinValueStrs[VarIdx2]))
					MinValueStrs[VarIdx2] = Neighbor[VarIdx2];
				if (MaxValueStrs[VarIdx2].empty() || Value > StrToFloat(MinValueStrs[VarIdx2]))
					MaxValueStrs[VarIdx2] = Neighbor[VarIdx2];
				}
			double dy = m_Best_y - y;
			xv2xss(Neighbor, xss);
			fprintf(f, "  ");
			fprintf(f, "%8.2g", dy);
			fprintf(f, "  %s", xss.c_str());
			fprintf(f, "\n");
			}
		}
	fprintf(f, "\n");
	fprintf(f, "%10.10s", "Min");
	fprintf(f, "  %10.10s", "Max");
	fprintf(f, "  %s", "Var");
	fprintf(f, "\n");
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		fprintf(f, "%10.4f", StrToFloat(MinValueStrs[VarIdx].c_str()));
		fprintf(f, "  %10.4f", StrToFloat(MaxValueStrs[VarIdx].c_str()));
		fprintf(f, "  %s", GetVarName(VarIdx));
		fprintf(f, "\n");
		}
	}

void Peaker::RunLatinClimb1()
	{
	string GlobalSpec;
	GetGlobalSpec(m_SpecLines, GlobalSpec);

	uint LatinBinCount = Peaker::SpecGetInt(GlobalSpec, "latin", UINT_MAX);
	asserta(LatinBinCount != UINT_MAX);
	RunLatin(LatinBinCount);
	HJ_RunHookeJeeves();
	WriteFinalResults(g_fLog);
	}

void Peaker::LoadTSV(const string &fn)
	{
	const uint nvar = GetVarCount();
	FILE *f = OpenStdioFile(fn);
	string line;
	vector<string> flds;
	// 1=y 2=why 3=dy 4=desc 5=varstr
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() == 5);
		double y = StrToFloat(flds[0]);
		const string &why = flds[1];
		const string &varstr = flds[4];
		vector<string> xv;
		xss2xv(varstr, xv);
		AppendResult(xv, y, why);
		}
	CloseStdioFile(f);
	}

static double get_random_float(double minx, double maxx)
	{
	const uint M = 1000;
	uint i = randu32()%M;
	double r = minx + (i*(maxx - minx)/(M-1));
	asserta(r >= minx && r <= maxx);
	return r;
	}

double Peaker::Bracket_FindNeighbor(
	uint VarIdx,
	bool Plus,
	double mindy_pct,
	double maxdy_pct,
	uint max_iters,
	string &value_str,
	double &rate)
	{
	vector<string> &saved_bestxv = m_Best_xv;
	const double saved_besty = m_Best_y;
	const string &start_str = saved_bestxv[VarIdx];
	double start_value = VarStrToFloat(VarIdx, start_str);

	if (rate == FLT_MAX) rate = 0.05;
	double upper_rate = FLT_MAX;
	double lower_rate = FLT_MAX;
	for (uint iter = 0; iter < max_iters; ++iter)
		{
		double old_rate = rate;
		vector<string> xv = saved_bestxv;

		double try_value = Plus ?
			start_value*(1 + rate) :
			start_value*(1 - rate);

		string try_str;
		Ps(try_str, "%.5g", try_value);
		xv[VarIdx] = try_str;
	
		double y = Evaluate(xv, "bktnbr");
		double dy = y - saved_besty;
		double pct = GetPct(abs(dy), y);
		ProgressLog("Bracket_FindNeighbor %c%s", pom(Plus), GetVarName(VarIdx));
		ProgressLog(" [%.5g]  y=%.5g", saved_besty, y);
		ProgressLog(" dy=%.5g (%.2f%%)", dy, pct);
		if (pct < mindy_pct)
			{
			lower_rate = rate;
			rate *= (1 + get_random_float(0.4, 0.8));
			if (upper_rate != FLT_MAX && rate > upper_rate)
				rate = (lower_rate + upper_rate)/2;
			ProgressLog(" iter %u/%u too small (min=%.3g%%) rate %.5g -> %.5g\n",
				iter+1, max_iters, mindy_pct, old_rate, rate);
			}
		else if (pct > maxdy_pct)
			{
			upper_rate = rate;
			rate /= (1 + get_random_float(0.4, 0.8));
			if (lower_rate != FLT_MAX && rate < lower_rate)
				rate = (lower_rate + upper_rate)/2;
			ProgressLog(" iter %u/%u too large (max=%.3g%%) rate %.5g -> %.5g\n",
				iter+1, max_iters, maxdy_pct, old_rate, rate);
			}
		else
			{
			ProgressLog(" OK\n");
			value_str = try_str;
			return y;
			}
		if (iter+1 == max_iters && y != saved_besty)
			{
			value_str = try_str;
			ProgressLog(" quit (iters)\n");
			return y;
			}
		}
	value_str = "";
	ProgressLog(" quit (dy=0)\n");
	return DBL_MAX;
	}

void Peaker::Bracket_FindNeighbors(
	double mindy_pct, double maxdy_pct, uint max_iters,
	vector<string> &plus_value_strs,
	vector<string> &minus_value_strs,
	vector<double> &plus_ys,
	vector<double> &minus_ys,
	vector<double> &rates)
	{
	asserta(mindy_pct < maxdy_pct);
	plus_value_strs.clear();
	minus_value_strs.clear();
	plus_ys.clear();
	minus_ys.clear();

	const uint VarCount = GetVarCount();

	plus_value_strs.resize(VarCount);
	minus_value_strs.resize(VarCount);
	plus_ys.resize(VarCount);
	minus_ys.resize(VarCount);

	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		ProgressLog("\n=== Neighbor %u/%u %s ===\n",
			VarIdx+1, VarCount, GetVarName(VarIdx));
		plus_ys[VarIdx] = 
			Bracket_FindNeighbor(
				VarIdx, true, mindy_pct, maxdy_pct,
				max_iters, plus_value_strs[VarIdx], rates[VarIdx]);

		minus_ys[VarIdx] = 
			Bracket_FindNeighbor(
				VarIdx, false, mindy_pct, maxdy_pct,
				max_iters, minus_value_strs[VarIdx], rates[VarIdx]);
		}
	}
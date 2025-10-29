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

double Peaker::rr(double lo, double hi) const
	{
	const uint M = 18409199;
	double r = (randu32()%(M+1))/double(M);
	double x = lo + r*(hi - lo);
	assert(x >= lo && x <= hi);
	return x;
	}

void VarSpec::Init(const vector<string> &Names,
	const vector<string> &Values,
	uint DefaultSigFig,
	bool RequireInitialValue)
	{
	m_SigFig = DefaultSigFig;
	const uint N = SIZE(Names);
	asserta(SIZE(Values) == N);
	for (uint i = 0; i < N; ++i)
		{
		const string &Name = Names[i];
		const string &Value = Values[i];

		if (Name == "var")
			m_Name = Value;
		else if (Name == "init")
			m_InitialValue = StrToFloat(Value);
		else if (Name == "min")
			m_Min = StrToFloat(Value);
		else if (Name == "max")
			m_Max = StrToFloat(Value);
		else if (Name == "delta")
			m_InitialDelta = StrToFloat(Value);
		else if (Name == "mindelta")
			m_MinDelta = StrToFloat(Value);
		else if (Name == "sigfig")
			m_SigFig = StrToUint(Value);
		else if (Name == "constant")
			{
			if (Value == "yes")
				m_Constant = true;
			else if (Value == "no")
				m_Constant = false;
			else
				Die("Invalid constant=%s", Value.c_str());
			}
		else
			Die("var %s spec name '%s'",
			  m_Name.c_str(), Name.c_str());
		}
	if (RequireInitialValue)
		asserta(m_InitialValue != DBL_MAX);
	if (m_Min == DBL_MAX)
		m_Min = m_InitialValue;
	if (m_Max == DBL_MAX)
		m_Max = m_InitialValue;
	if (m_InitialDelta == DBL_MAX)
		m_InitialDelta = fabs(m_InitialValue/10);
	if (m_MinDelta == DBL_MAX)
		m_MinDelta = m_InitialDelta/10;
	}

const char *Peaker::GetVarName(uint VarIdx) const
	{
	const VarSpec &Spec = GetVarSpec(VarIdx);
	return Spec.m_Name.c_str();
	}

void Peaker::Round_xs(const vector<double> &xv, vector<double> &Rounded_xv) const
	{
	Rounded_xv.clear();
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv) == VarCount);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		double Value = xv[VarIdx];
		uint SigFig = GetSigFig(VarIdx);
		double RoundedValue = GetRounded(Value, SigFig);
		Rounded_xv.push_back(RoundedValue);
		}
	}

bool Peaker::Cmp_xs(const vector<double> &xs1, const vector<double> &xs2) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xs1) == VarCount);
	asserta(SIZE(xs2) == VarCount);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const VarSpec &Spec = GetVarSpec(VarIdx);
		double x1 = xs1[VarIdx];
		double x2 = xs2[VarIdx];
		if (!feq(x1, x2))
			return false;
		}
	return true;
	}

uint Peaker::Find_xs(const vector<double> &xv) const
	{
	for (uint xIdx = 0; xIdx < SIZE(m_xvs); ++xIdx)
		{
		const vector<double> &xs2 = m_xvs[xIdx];
		if (Cmp_xs(xv, xs2))
			return xIdx;
		}
	return UINT_MAX;
	}

void Peaker::LogLatinBins() const
	{
	uint VarCount = GetVarCount();
	Log("\n");
	Log("LogLatinBins [%u]\n", m_LatinBinCount);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const VarSpec &Spec = GetVarSpec(VarIdx);
		const char *Name = GetVarName(VarIdx);
		if (Spec.m_Constant)
			{
			Log("%s constant %.3g\n", Name, Spec.m_InitialValue);
			continue;
			}
		Log("%s %.3g .. %.3g\n", Name, Spec.m_Min, Spec.m_Max);
		Log("  ");
		for (uint BinIdx = 0; BinIdx < m_LatinBinCount; ++BinIdx)
			Log(" %.3g", GetLatinValueByBinIdx(VarIdx, BinIdx, m_LatinBinCount));
		Log("\n");
		}
	}

double Peaker::GetLatinBinWidth(uint VarIdx) const
	{
	const VarSpec &Spec = GetVarSpec(VarIdx);
	if (Spec.m_Constant)
		return 0;
	asserta(Spec.m_Min != DBL_MAX);
	asserta(Spec.m_Max != DBL_MAX);
	asserta(Spec.m_Min < Spec.m_Max);
	double BinWidth = (Spec.m_Max - Spec.m_Min)/m_LatinBinCount;
	return BinWidth;
	}

double Peaker::GetLatinValueByBinIdx(uint VarIdx, uint BinIdx, uint BinCount) const
	{
	asserta(BinCount > 0);
	const VarSpec &Spec = GetVarSpec(VarIdx);
	if (Spec.m_Constant)
		return Spec.m_InitialValue;
	asserta(Spec.m_Min != DBL_MAX);
	asserta(Spec.m_Max != DBL_MAX);
	asserta(Spec.m_Min < Spec.m_Max);
	double BinWidth = (Spec.m_Max - Spec.m_Min)/BinCount;
	asserta(feq(BinWidth*BinCount, Spec.m_Max - Spec.m_Min));
	double BinLo = Spec.m_Min + BinWidth*BinIdx;
	const uint M = 3141592;
	double r = double(randu32()%M)/(M-1);
	asserta(r >= 0 && r <= 1);
	double Value = BinLo + r*BinWidth;
	return Value;
	}

void Peaker::GetLatinHypercube(vector<vector<double> > &xvs)
	{
	xvs.clear();
	if (m_LatinBinCount == 0)
		return;

	const uint VarCount = GetVarCount();
	vector<vector<uint> > IdxMx;
	GetLatinHypercubeIdxs(VarCount, m_LatinBinCount, IdxMx);
	asserta(SIZE(IdxMx) == m_LatinBinCount);

	xvs.resize(m_LatinBinCount);
	for (uint k = 0; k < m_LatinBinCount; ++k)
		{
		vector<double> &xv = xvs[k];
		for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
			{
			uint BinIdx = IdxMx[k][VarIdx];
			double Value = GetLatinValueByBinIdx(VarIdx, BinIdx, m_LatinBinCount);
			xv.push_back(Value);
			}
		}
	}

const VarSpec &Peaker::GetVarSpec(uint VarIdx) const
	{
	asserta(VarIdx < SIZE(m_VarSpecs));
	return *m_VarSpecs[VarIdx];
	}

bool Peaker::VarIsConstant(uint VarIdx) const
	{
	return GetVarSpec(VarIdx).m_Constant;
	}

double Peaker::GetMinDelta(uint VarIdx) const
	{
	return GetVarSpec(VarIdx).m_MinDelta;
	}

void Peaker::InitDeltas()
	{
	m_Deltas.clear();
	const uint VarCount = GetVarCount();
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const VarSpec &Spec = GetVarSpec(VarIdx);
		double Delta = Spec.m_InitialDelta;
		m_Deltas.push_back(Delta);
		}
	}

void Peaker::LogDeltas() const
	{
	const uint VarCount = GetVarCount();

	Log("Deltas: ");
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const VarSpec &Spec = GetVarSpec(VarIdx);
		double InitialDelta = Spec.m_InitialDelta;
		Log(" %s(%.3g,%.3g)", GetVarName(VarIdx), InitialDelta, m_Deltas[VarIdx]);
		}
	Log("\n");
	}

uint Peaker::GetSigFig(uint VarIdx) const
	{
	const VarSpec &Spec = GetVarSpec(VarIdx);
	if (Spec.m_SigFig != UINT_MAX)
		return Spec.m_SigFig;
	return m_SigFig;
	}

const char *Peaker::VarsToStr(const vector<double> &xv, string &s,
							  const string sep) const
	{
	const uint VarCount = GetVarCount();
	asserta(SIZE(xv) == VarCount);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		double x = xv[VarIdx];
		asserta(x != DBL_MAX);
		asserta(!isnan(x));
		string Tmp;
		VarToStr(x, VarIdx, Tmp);
		if (VarIdx > 0)
			s += sep;
		Psa(s, "%s=%s", GetVarName(VarIdx), Tmp.c_str());
		}
	return s.c_str();
	}

const char *Peaker::VarToStr(double x, uint VarIdx, string &s) const
	{
	asserta(x != DBL_MAX);
	asserta(!isnan(x));
	uint SigFig = GetSigFig(VarIdx);

	string Fmt;
	Ps(Fmt, "%%.%ug", SigFig);
	Ps(s, Fmt.c_str(), x);
	return s.c_str();
	}

double Peaker::Calc(const vector<double> &xv)
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
	if (m_ptrWhoAmI != 0)
		WhoAmI = *m_ptrWhoAmI;
	s = ParentStr + "/" + WhoAmI;
	}

void Peaker::AppendChildResult(const vector<double> &xv, double y, 
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
		string VarStr;
		VarsToStr(xv, VarStr);
		fprintf(m_fTsv, "%.6g\t%s%s\t%s\n",
			y, improved.c_str(), why.c_str(), VarStr.c_str());
		fflush(m_fTsv);
		}
	}

double Peaker::Evaluate(const vector<double> &axv, const string &why)
	{
	const double Saved_Best_y = m_Best_y;
	vector<double> xv;
	Round_xs(axv, xv);
	uint Idx = Find_xs(xv);
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
	string VarsStr;
	VarsToStr(xv, VarsStr, ";");
	ProgressLog("[%.6g](%+.3g) %s %s %s\n",
		m_Best_y, dy, Path.c_str(), why.c_str(), VarsStr.c_str());
	string improved;
	if (dy > 0)
		improved = ">>";
	if (m_fTsv != 0 && y != DBL_MAX)
		{
		fprintf(m_fTsv, "%.6g\t%s%s:%s\t%s\n",
			y, improved.c_str(), Path.c_str(), why.c_str(), VarsStr.c_str());
		fflush(m_fTsv);
		}
	return y;
	}

void Peaker::Init(const vector<string> &SpecLines, PTR_EVAL_FUNC EF)
	{
	Clear();
	m_EvalFunc = EF;
	vector<string> Fields;
	const uint LineCount = SIZE(SpecLines);
	m_SkipInit = false;
	for (uint LineIdx = 0; LineIdx < LineCount; ++LineIdx)
		{
		const string &Line = SpecLines[LineIdx];
		if (SIZE(Line) == 0 || Line[0] == '#')
			continue;
		Split(Line, Fields, '\t');
		const uint n = SIZE(Fields);
		if (n == 0)
			continue;
		vector<string> Names;
		vector<string> Values;
		for (uint i = 0; i < n; ++i)
			{
			const string &Field = Fields[i];
			vector<string> Fields2;
			Split(Field, Fields2, '=');
			if (SIZE(Fields2) != 2)
				Die("Expected name=value '%s'", Field.c_str());
			const string &Name =Fields2[0];
			const string &Value = Fields2[1];
			Names.push_back(Name);
			Values.push_back(Value);
			}
		if (SIZE(Names) == 0)
			continue;
		if (Names[0] == "var")
			{
			VarSpec *Spec = new VarSpec;
			Spec->Init(Names, Values, m_SigFig, !m_SkipInit);
			m_VarSpecs.push_back(Spec);
			}
		else
			{
			asserta(SIZE(Names) == 1);
			asserta(SIZE(Values) == 1);
			const string &Name = Names[0];
			const string &Value = Values[0];
			if (Name == "dy")
				m_Target_dy = StrToFloat(Value);
			else if (Name == "cmd")
				m_Cmd = Value;
			else if (Name == "mindy")
				m_Min_dy = StrToFloat(Value);
			else if (Name == "maxdy")
				m_Max_dy = StrToFloat(Value);
			else if (Name == "sigfig")
				m_SigFig = StrToUint(Value);
			else if (Name == "latin")
				m_LatinBinCount = StrToUint(Value);
			else if (Name == "hj_iters")
				m_HJ_MaxIters = StrToUint(Value);
			else if (Name == "minh")
				m_Min_Height = StrToFloat(Value);
			else if (Name == "qdir")
				m_QueueDir = Value;
			else if (Name == "skipinit")
				{
				if (Value == "yes")
					m_SkipInit = true;
				else if (Value == "no")
					m_SkipInit = false;
				else
					Die("Invalid skipinit=%s", Value.c_str());
				}
			else
				Die("Invalid param name '%s'", Name.c_str());
			}
		}
	}

void Peaker::LogState() const
	{
	Die("TODO");
	}

void Peaker::GetBestVars(vector<double> &xs) const
	{
	xs = m_Best_xv;
	}

void Peaker::GetBestVarStr(string &s) const
	{
	VarsToStr(m_Best_xv, s);
	}

double Peaker::GetRounded(double x, uint SigFig) const
	{
	string Fmt;
	Ps(Fmt, "%%.%ug", SigFig);
	string Tmp;
	Ps(Tmp, Fmt.c_str(), x);
	double Rounded_x = atof(Tmp.c_str());
	return Rounded_x;
	}

void Peaker::RunInitialValues()
	{
	const uint VarCount = GetVarCount();
	vector<double> xv;
	uint n = 0;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const VarSpec &Spec = GetVarSpec(VarIdx);
		double InitialValue = Spec.m_InitialValue;
		if (InitialValue != DBL_MAX)
			{
			xv.push_back(InitialValue);
			++n;
			}
		}
	if (n == 0)
		return;
	if (n < VarCount)
		Die("Missing %u / %u initial values",
		  VarCount - n, VarCount);
	Evaluate(xv, "init");
	}

void Peaker::RunLatin(vector<vector<double> > &xvs, vector<double> &ys)
	{
	xvs.clear();
	ys.clear();
	if (m_LatinBinCount == 0)
		return;

	const uint VarCount = GetVarCount();
	GetLatinHypercube(xvs);
	const uint n = SIZE(xvs);
	asserta(n > 0);
#if	1
	{
	const uint m = GetVarCount();
	Log("Latin=%u\n", n);
	for (uint i = 0; i < n; ++i)
		{
		Log("[%5u] ", i);
		const vector<double> &xv = xvs[i];
		asserta(SIZE(xv) == m);
		for (uint j = 0; j < m; ++j)
			Log(" %5.5s=%8.3g", GetVarName(j), xv[j]);
		Log("\n");
		}
	}
#endif
	for (uint i = 0; i < n; ++i)
		{
		string why;
		Ps(why, "latin %u/%u", i+1, n);
		double y = Evaluate(xvs[i], why);
		ys.push_back(y);
		}
	asserta(SIZE(xvs) == SIZE(ys));
	}

void Peaker::CleanQueue()
	{
	asserta(m_QueueDir != "");
	asserta(m_QueueDir != "/");
	asserta(m_QueueDir != ".");
	vector<string> FNs;
	vector<bool> IsSubDirs;
	mylistdir(m_QueueDir, FNs, IsSubDirs);
	const uint n = SIZE(FNs);
	for (uint i = 0; i < n; ++i)
		{
		if (IsSubDirs[i])
			continue;
		const string &FN = FNs[i];
		if (EndsWith(FN, ".request") ||
			EndsWith(FN, ".request.x") ||
			EndsWith(FN, ".tmp") ||
			EndsWith(FN, ".result") ||
			EndsWith(FN, ".done"))
			{
			Progress("delete %s\n", FN.c_str());
			DeleteStdioFile(m_QueueDir + FN);
			}
		}
	}

void Peaker::Run()
	{
	RunInitialValues();
	vector<vector<double> > latin_xvs;
	vector<double> ys;
	RunLatin(latin_xvs, ys);
	HJ_RunHookeJeeves();
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

void Peaker::LogSpecs() const
	{
	const uint VarCount = GetVarCount();
	Log("%u vars\n", VarCount);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const VarSpec &Spec = GetVarSpec(VarIdx);
		Log("var=%s", Spec.m_Name.c_str());
		Log("\tinit=%.3g", Spec.m_InitialValue);
		if (Spec.m_Constant)
			Log("\tconstant=yes");
		else
			{
			Log("\tmin=%.3g", Spec.m_Min);
			Log("\tmax=%.3g", Spec.m_Max);
			Log("\tsigfig=%u", GetSigFig(VarIdx));
			}
		Log("\n");
		}
	}

void Peaker::RunNestedLatin(uint TopN)
	{
	ProgressLog("RunNestedLatin(%u)\n", TopN);
	vector<vector<double> > latin_xvs;
	vector<double> ys;
	RunLatin(latin_xvs, ys);

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
		const vector<double> &xv = m_xvs[i];
		double yseed = ys[i];
		asserta(SIZE(xv) == VarCount);
		Log("RUN_NESTED_LATIN %u/%u y=%.3g\n", k+1, N, ys[i]);

		string nest;
		Ps(nest, "NestedLatin%u", k+1);
		Peaker P(this, &nest);
		P.m_LatinBinCount = m_LatinBinCount;
		P.m_Target_dy = m_Target_dy;
		P.m_Min_dy = m_Min_dy;
		P.m_Max_dy = m_Max_dy;
		P.m_Min_Height = m_Min_Height;
		P.m_SigFig = m_SigFig;
		P.m_EvalFunc = m_EvalFunc;
		P.m_fTsv = m_fTsv;
		for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
			{
			double Value = xv[VarIdx];
			double BinWidth = GetLatinBinWidth(VarIdx);
			const VarSpec &Spec = GetVarSpec(VarIdx);
			VarSpec &SubSpec = *new VarSpec;
			SubSpec = Spec;
			SubSpec.m_InitialValue = Value;
			SubSpec.m_Min = Value - BinWidth/2;
			if (SubSpec.m_Min < 0)
				SubSpec.m_Min = 0;
			SubSpec.m_Max = Value + BinWidth/2;
			P.m_VarSpecs.push_back(&SubSpec);
			}

		Log("\n");
		Log("NESTED_LATIN\n");
		P.LogSpecs();
		vector<vector<double> > Child_xvs;
		vector<double> Child_ys;
		P.RunLatin(Child_xvs, Child_ys);
		const uint J = SIZE(Child_xvs);
		asserta(SIZE(P.m_whys) == J);
		asserta(SIZE(Child_ys) == J);
		Log("NESTED_LATIN_DONE %u/%u yseed0 %.3g yseed=%.3g best=%.3g J=%u\n",
			k+1, N, yseed0, yseed, P.m_Best_y, J);
		for (uint j = 0; j < J; ++j)
			{
			const vector<double> &Child_xv = Child_xvs[j];
			double y = Child_ys[j];
			const string &Child_why = P.m_whys[j];
			AppendChildResult(Child_xv, y, nest, Child_why);
			}
		}
	}

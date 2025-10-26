#include "myutils.h"
#include "peaker.h"
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
  const vector<string> &Values)
	{
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
	asserta(m_InitialValue != DBL_MAX);
	if (m_Min == DBL_MAX)
		m_Min = m_InitialValue;
	if (m_Max == DBL_MAX)
		m_Max = m_InitialValue;
	if (m_InitialDelta == DBL_MAX)
		m_InitialDelta = fabs(m_InitialValue/10);
	if (m_MinDelta == DBL_MAX)
		m_MinDelta = m_InitialDelta/10;
	if (m_SigFig == UINT_MAX)
		m_SigFig = 3;
	}

const char *Peaker::GetVarName(uint VarIdx) const
	{
	const VarSpec &Spec = GetVarSpec(VarIdx);
	return Spec.m_Name.c_str();
	}

// xv are equivalent if all diffs are <= VarSpec.MinDelta
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

double Peaker::Get_y(uint xIdx) const
	{
	const uint N = SIZE(m_ys);
	asserta(N == SIZE(m_xvs));
	asserta(xIdx < N);
	return m_ys[xIdx];
	}

const vector<double> &Peaker::Get_xv(uint xIdx) const
	{
	asserta(xIdx < SIZE(m_xvs));
	return m_xvs[xIdx];
	}

double Peaker::Get_x(uint xIdx, uint VarIdx) const
	{
	const vector<double> &xv = Get_xv(xIdx);
	asserta(VarIdx < SIZE(xv));
	double x = xv[VarIdx];
	return x;
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

uint Peaker::Add_xv(const vector<double> &xv)
	{
	uint xIdx = Find_xs(xv);
	if (xIdx != UINT_MAX)
		return xIdx;
	xIdx = SIZE(m_xvs);

	const uint VarCount = GetVarCount();
	vector<double> Rounded_xv;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		double x = xv[VarIdx];
		uint SigFig = GetSigFig(VarIdx);
		double Rounded_x = GetRounded(x, SigFig);
		Rounded_xv.push_back(Rounded_x);
		}

	m_xvs.push_back(Rounded_xv);
	m_ys.push_back(DBL_MAX);
	return xIdx;
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
	double BinWidth = (Spec.m_Max - Spec.m_Min)/(BinCount - 1);
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
	for (uint BinIdx = 0; BinIdx < m_LatinBinCount; ++BinIdx)
		{
		vector<double> &xv = xvs[BinIdx];
		for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
			{
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
	if (m_Progress)
		Progress("Deltas: ");
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const VarSpec &Spec = GetVarSpec(VarIdx);
		double InitialDelta = Spec.m_InitialDelta;
		Log(" %s(%.3g,%.3g)", GetVarName(VarIdx), InitialDelta, m_Deltas[VarIdx]);
		if (m_Progress)
			Progress(" %s(%.3g,%.3g)", GetVarName(VarIdx), InitialDelta, m_Deltas[VarIdx]);
		}
	Log("\n");
	if (m_Progress)
		Progress("\n");
	}

uint Peaker::GetSigFig(uint VarIdx) const
	{
	const VarSpec &Spec = GetVarSpec(VarIdx);
	if (Spec.m_SigFig != UINT_MAX)
		return Spec.m_SigFig;
	if (m_SigFig != UINT_MAX)
		return m_SigFig;
	return m_DefaultSigFig;
	}

const char *Peaker::VarsToStr(const vector<double> &xv, string &s,
							  const string sep) const
	{
	const uint VarCount = GetVarCount();
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

double Peaker::Evaluate(uint xIdx, bool UnsetOk)
	{
	if (xIdx == UINT_MAX)
		{
		if (UnsetOk)
			return DBL_MAX;
		Die("Evaluate(*)");
		}
	asserta(xIdx < SIZE(m_ys));
	asserta(xIdx < SIZE(m_xvs));
	const vector<double> &xv = m_xvs[xIdx];
	string VarsStr;
	VarsToStr(xv, VarsStr);
	if (m_ys[xIdx] != DBL_MAX)
		{
		double y = m_ys[xIdx];
		++m_EvaluateCacheHits;
		return y;
		}

	double y = Calc(xv);
	if (y == DBL_MAX)
		return DBL_MAX;
	m_ys[xIdx] = y;
	double dy = DBL_MAX;
	if (m_Best_y == DBL_MAX)
		dy = fabs(y);
	else
		dy = y - m_Best_y;
	Log("y=%.6g/%+.3g Evaluate(%s)\n", y, dy, VarsStr.c_str());
	if (dy > 0)
		{
		m_Best_xIdx = xIdx;
		m_Best_y = y;
		m_Best_ys.push_back(y);
		m_Best_EvalIdxs.push_back(m_EvalIdx);
		m_Best_xIdxs.push_back(xIdx);
		m_LastImprovedEvalIdx = m_EvalIdx;
		ProgressLogSummary();
		}
	++m_EvalIdx;
	if (m_fTsv != 0)
		{
		fprintf(m_fTsv, "%.6g\t%s\n", y, VarsStr.c_str());
		fflush(m_fTsv);
		}
	return y;
	}

void Peaker::LogPair(uint xIdx1, uint xIdx2) const
	{
	Log("Pair: ");
	const vector<double> &xv1 = Get_xv(xIdx1);
	const vector<double> &xv2 = Get_xv(xIdx2);
	const uint VarCount = GetVarCount();
	uint ChangeCount = 0;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		double x1 = xv1[VarIdx];
		double x2 = xv2[VarIdx];
		if (x1 != x2)
			{
			++ChangeCount;
			Log(" %s %.6g(%+.3g)", GetVarName(VarIdx), x1, x2 - x1);
			}
		}
	Log(" %u changes\n", ChangeCount);
	}

void Peaker::Init(const vector<string> &SpecLines,
  PTR_EVAL_FUNC EF)
	{
	Clear();
	m_EvalFunc = EF;
	vector<string> Fields;
	const uint LineCount = SIZE(SpecLines);
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
			Spec->Init(Names, Values);
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
			else
				Die("Invalid param name '%s'", Name.c_str());
			}
		}
	}

void Peaker::LogState() const
	{
	const uint VarCount = GetVarCount();
	vector<double> xv;
	if (m_Best_xIdx != UINT_MAX)
		xv = Get_xv(m_Best_xIdx);

	Log("\n");
	Log("#%u", m_EvalIdx);
	Log(", last+ %u", m_LastImprovedEvalIdx);
	Log(", cache %u", m_EvaluateCacheHits);
	Log(", target dy=%.3g", m_Target_dy);
	Log("\n");
	Log("       Var       Delta        Best\n");
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		Log("%10.10s", GetVarName(VarIdx));
		Log("  %10.10s", FloatToStr(m_Deltas[VarIdx]));
		if (m_Best_xIdx != UINT_MAX)
			Log("  %10.6g", xv[VarIdx]);
		Log("\n");
		}
	}

void Peaker::ProgressLogSummary() const
	{
	string TmpStr;
	Log("%s ", GetElapsedTimeStr(TmpStr));
	if (m_Progress)
		Progress("%s ", GetElapsedTimeStr(TmpStr));
	if (m_Best_xIdx == UINT_MAX)
		{
		Log(" (no evals)\n");
		if (m_Progress)
			Progress(" (no evals)\n");
		return;
		}

	if (m_LastImprovedEvalIdx == m_EvalIdx)
		{
		Log(">>", m_EvalIdx);
		if (m_Progress)
			Progress(">>", m_EvalIdx);
		}
	else
		{
		Log("#%u", m_EvalIdx);
		Log(" ~%u", m_EvalIdx - m_LastImprovedEvalIdx);
		if (m_Progress)
			{
			Progress("#%u", m_EvalIdx);
			Progress(" ~%u", m_EvalIdx - m_LastImprovedEvalIdx);
			}
		}
	if (m_Msg != "")
		{
		Log(" %s", m_Msg.c_str());
		if (m_Progress)
			Progress(" %s", m_Msg.c_str());
		}

	string s;
	GetBestVarStr(s);
	Log(" [%.6g]", m_Best_y);
	Log(" %s", s.c_str());
	Log("\n");
	if (m_Progress)
		{
		Progress(" [%.6g]", m_Best_y);
		Progress(" %s", s.c_str());
		Progress("\n");
		}
	}

void Peaker::GetBestVarStr(string &s) const
	{
	if (m_Best_xIdx == UINT_MAX)
		{
		s = "(none)";
		return;
		}
	const vector<double> &xv = Get_xv(m_Best_xIdx);
	VarsToStr(xv, s);
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
	uint TmpIdx = Add_xv(xv);
	m_Msg = "Inits";
	Evaluate(TmpIdx);
	}

void Peaker::RunLatin()
	{
	if (m_LatinBinCount == 0)
		return;

	const uint VarCount = GetVarCount();
	vector<vector<double> > xvs;
	GetLatinHypercube(xvs);
	const uint n = SIZE(xvs);
	asserta(n > 0);
#if	0
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
		Ps(m_Msg, "Latin%u/%u", i+1, n);
		const vector<double> &xv = xvs[i];
		uint uIdx = Add_xv(xv);
		Log("Latin %u/%u\n", i+1, n);
		Evaluate(uIdx);
		}
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

const vector<double> &Peaker::GetBestParams() const
	{
	asserta(m_Best_xIdx < SIZE(m_xvs));
	return m_xvs[m_Best_xIdx];
	}

void Peaker::LatinLoop()
	{
	for (;;)
		RunLatin();
	}

void Peaker::Run()
	{
	RunInitialValues();
	RunLatin();
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

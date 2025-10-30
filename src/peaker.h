#pragma once

typedef double (*PTR_EVAL_FUNC)(const vector<string> &xv);

static const uint MIN_RATE = 1;
static const uint MED_RATE = 3;
static const uint MAX_RATE = 5;

class Peaker
	{
public:
	static FILE *m_fTsv;

public:
	string m_SpecFN;
	vector<string> m_SpecLines;

	const Peaker *m_Parent = 0;
	string m_WhoAmI;
	string m_GlobalSpec;
	vector<string> m_VarSpecs;
	vector<string> m_VarNames;
	uint m_EvaluateCacheHits = 0;
	uint m_RequestIdx = 0;
	PTR_EVAL_FUNC m_EvalFunc = 0;

	vector<vector<string> > m_xvs;		// Values tried
	vector<double> m_ys;				// Results
	vector<string> m_whys;				// Reason for Eval
	double m_Best_y = DBL_MAX;
	vector<string> m_Best_xv;
	time_t m_LastImprovedTime = 0;
	uint m_LastImprovedEvalCount = 0;

	vector<uint> m_Rates;

// Hooke-Jeeves parameters
	uint m_HJ_MaxExtendIters = 100;
	uint m_HJ_MaxIters = 100;

// Hooke-Jeeves state
	uint m_HJ_Direction = UINT_MAX;		// current axis
	bool m_HJ_ExtendPlus = false;	// +/- direction

public:
	Peaker(const Peaker *Parent, const string &WhoAmI)
		{
		m_Parent = Parent;
		m_WhoAmI = WhoAmI;
		}

	void Clear()
		{
		Die("Peaker::Clear()");
		m_Parent = 0;
		m_WhoAmI.clear();
		m_SpecFN.clear();
		m_SpecLines.clear();
		m_VarNames.clear();
		m_VarSpecs.clear();
		m_EvaluateCacheHits = 0;
		m_RequestIdx = 0;
		m_EvalFunc = 0;

		m_xvs.clear();
		m_ys.clear();
		m_whys.clear();
		m_Best_y = DBL_MAX;
		m_LastImprovedTime = 0;
		m_LastImprovedEvalCount = 0;
		m_Best_xv.clear();
		m_Rates.clear();
		}

public:
	void Init(const vector<string> &SpecLines, PTR_EVAL_FUNC EF);
	uint Find_xv(const vector<string> &xv) const;
	void WriteStatusPage(FILE *f) const;
	void LogSpec() const;
	void InitRates();
	const char *GetVarName(uint VarIdx) const;
	uint GetVarCount() const { return SIZE(m_VarSpecs); }
	void GetLatinHypercube(vector<vector<string> > &xvs);
	double Evaluate(const vector<string> &xv, const string &why);
	double Calc(const vector<string> &xv);
	void Run();
	void LogState() const;
	void xv2str(const vector<string> &xv, string &xstr) const;
	void GetBestVarStr(string &s) const;
	void GetBestVars(vector<double> &xs) const;
	void RunInitialValues();
	void RunLatin();
	double rr(double lo, double hi) const;
	bool VarIsConstant(uint VarIdx) const;
	double GetBinWidth(uint VarIdx) const;
	void NormalizeWeights(const vector<string> &xv,
		vector<string> &Normalized_xv) const;

	const string &GetVarSpec(uint VarIdx) const;

	void VarFloatToStr(uint VarIdx, double Value, string &s) const;
	double VarStrToFloat(uint VarIdx, const string &ValueStr) const;

	double VarSpecGetFloat(uint VarIdx, const string &Name, double Default) const;
	uint VarSpecGetInt(uint VarIdx, const string &Name, uint Default) const;
	bool VarSpecGetBool(uint VarIdx, const string &Name, bool Default) const;
	void VarSpecGetStr(uint VarIdx, const string &Name, 
		string &Str, const string &Default) const;

	void GetGlobalStr(const string &Name, string &s, const string &Default) const;
	double GetGlobalFloat(const string &Name, double Default) const;
	uint GetGlobalInt(const string &Name, uint Default) const;
	bool GetGlobalBool(const string &Name, bool Default) const;

	double GetLatinValueByBinIdx(uint VarIdx, uint BinIdx, uint BinCount) const;
	void RunNestedLatin(uint TopN);
	void GetPeakerPathStr(string &s) const;
	void AppendChildResult(const vector<string> &xv, double y,
		const string &childname, const string &why);
	void IncreaseRate(uint VarIdx);
	void DecreaseRate(uint VarIdx);
	void NormalizeVarStr(uint VarIdx, const string &Str,
		string &NormalizedStr) const;

	void DeltaVar(uint VarIdx, bool Plus, const string &OldStr,
		string &NewStr);

// Hooke-Jeeves
	void HJ_RunHookeJeeves();
	void HJ_Explore();
	void HJ_Extend();
	double HJ_TryDelta(const string &reason,
		const vector<string> &Start_xv, uint VarIdx, bool Plus);

public:
	static double GetRoundedValue(double x, uint SigFig);
	static void GetRoundedStr(double x, uint SigFig, string &Str);
	static void GetLatinHypercubeIdxs(uint VarCount, uint BinCount,
		vector<vector<uint> > &IdxMx);
	static void SpecReplaceStr(string &Spec, const string &Name, const string &NewValue);
	static double SpecGetFloat(const string &Spec, const string &Name, double Default);
	static void SpecGetStr(const string &Spec, const string &Name,
		string &Str, const string &Default);
	static bool SpecGetBool(const string &Spec, const string &Name, bool Default);
	static uint SpecGetInt(const string &Spec, const string &Name, uint Default);
	static double GetRateFactor(uint Rate, bool Plus);
	static double GetIncreaseRateFactor(uint Rate);
	static double GetDecreaseRateFactor(uint Rate);
	static void IncFloat(const string &OldStr, bool Plus, string &NewStr);
	static void ParseEStr(const string &EStr, string &Mantissa, string &Exponent);
	static uint GetSigFig(const string &EStr);
	static bool AllNines(const string &s);
	static bool AllZeros(const string &s);
	static bool OneZeros(const string &s);
	};

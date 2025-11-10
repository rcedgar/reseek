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
	string m_InitParams;

	const Peaker *m_Parent = 0;
	string m_Name;
	string m_GlobalSpec;
	vector<string> m_VarSpecs;
	vector<string> m_VarNames;
	uint m_GlobalVarRateFactorIdx = UINT_MAX;
	uint m_EvaluateCacheHits = 0;
	uint m_RequestIdx = 0;
	PTR_EVAL_FUNC m_EvalFunc = 0;

	vector<vector<string> > m_xvs;		// Values tried
	vector<double> m_ys;				// Results
	vector<string> m_whys;				// Reason for Eval
	double m_Best_y = DBL_MAX;
	vector<string> m_Best_xv;
	string m_Best_desc;
	time_t m_LastImprovedTime = 0;
	uint m_LastImprovedEvalCount = 0;

	vector<string> m_ChildNames;
	vector<vector<string> > m_Best_xvs;
	vector<double> m_Best_ys;
	vector<string> m_Best_descs;

// Hooke-Jeeves parameters
	uint m_HJ_MaxExtendIters = 100;
	uint m_HJ_MaxIters = 100;

// Hooke-Jeeves state
	uint m_HJ_Direction = UINT_MAX;		// current axis
	bool m_HJ_ExtendPlus = false;	// +/- direction

public:
	Peaker(const Peaker *Parent, const string &Name)
		{
		m_Parent = Parent;
		m_Name = Name;
		}

public:
	void Init(const vector<string> &SpecLines, PTR_EVAL_FUNC EF);
	void WriteFinalResults(FILE *f) const;
	void WriteStatusPage(FILE *f) const;
	void LogSpec() const;
	void LogState() const;
	void InitRates();

	// Child
	void GetPeakerPathStr(string &s) const;
	Peaker *MakeChild(const string &Name) const;
	void AppendResult(const vector<string> &xv, double y, const string &why);
	void AppendChildResults(const Peaker &Child);

	uint GetVarIdx(const string &Name, bool FailOk = false) const;
	uint GetVarCount() const { return SIZE(m_VarSpecs); }

	double Evaluate(const vector<string> &xv, const string &why);
	double Calc(const vector<string> &xv);
	void RunLatinClimb1();
	void RunLatin(uint BinCount);
	void NormalizeWeights(const vector<string> &xv,
		vector<string> &Normalized_xv) const;
	void DeltaVar(uint VarIdx, bool Plus, const string &OldStr,
		string &NewStr);

	bool GetNearestNeighbor(const vector<string> &xv,
		uint VarIdx, bool Plus, vector<string> &Neighbor_xv, double &y) const;
	void WriteFinalPeak(FILE *f) const;

	// Latin
	void GetLatinHypercube(uint BinCount, vector<vector<string> > &xvs) const;
	double GetLatinValueByBinIdx(uint VarIdx, uint BinIdx, uint BinCount) const;

	// Global spec
	void GetGlobalStr(const string &Name, string &s, const string &Default) const;
	double GetGlobalFloat(const string &Name, double Default) const;
	uint GetGlobalInt(const string &Name, uint Default) const;
	bool GetGlobalBool(const string &Name, bool Default) const;

	// Per-var spec
	const string &GetVarSpec(uint VarIdx) const;
	double VarSpecGetFloat(uint VarIdx, const string &Name, double Default) const;
	uint VarSpecGetInt(uint VarIdx, const string &Name, uint Default) const;
	bool VarSpecGetBool(uint VarIdx, const string &Name, bool Default) const;
	void VarSpecGetStr(uint VarIdx, const string &Name,
		string &Str, const string &Default) const;

	const char *GetVarName(uint VarIdx) const;
	bool VarIsConstant(uint VarIdx) const;
	double GetBinWidth(uint VarIdx) const;

	// xv vector of strings, xss is semi-colon string
	void xv2values(const vector<string> &xv, vector<double> &Values) const;
	void xss2xv(const string &xstr, vector<string> &xv) const;
	const char *xv2xss(const vector<string> &xv, string &xstr) const;
	uint Find_xv(const vector<string> &xv) const;

	// Convert float <-> str for one var
	void VarFloatToStr(uint VarIdx, double Value, string &s) const;
	double VarStrToFloat(uint VarIdx, const string &ValueStr) const;

	void NormalizeVarStr(uint VarIdx, const string &Str,
		string &NormalizedStr) const;
	double GetEuclideanDist(const vector<string> &xv1,
		const vector<string> &xv2) const;
	double GetEuclideanDist(const vector<double> &xv1,
		const vector<double> &xv2) const;

	void GetTopEvalIdxs(const uint N, vector<uint> &Idxs) const;

// Hooke-Jeeves
	void HJ_RunHookeJeeves();
	void HJ_Explore();
	void HJ_Extend();
	bool HJ_Iter();
	double HJ_TryDelta(const string &reason,
		const vector<string> &Start_xv, uint VarIdx, bool Plus,
		vector<string> &Try_xv);

	bool ReduceGlobalRateFactor();
	double GetGlobalRateFactor();
	double GetRateFactor(bool Plus);
	double GetIncreaseRateFactor();
	double GetDecreaseRateFactor();

public:
	static void GetGlobalSpec(const vector<string> &SpecLines,
		string &GlobalSpec);
	static double GetRoundedValue(double x, uint SigFig);
	static void GetRoundedStr(double x, uint SigFig, string &Str);
	static void GetLatinHypercubeIdxs(uint VarCount, uint BinCount,
		vector<vector<uint> > &IdxMx);
	static double SpecGetFloat(const string &Spec, const string &Name, double Default);
	static void SpecGetStr(const string &Spec, const string &Name,
		string &Str, const string &Default);
	static bool SpecGetBool(const string &Spec, const string &Name, bool Default);
	static uint SpecGetInt(const string &Spec, const string &Name, uint Default);
	static void IncFloat(const string &OldStr, bool Plus, string &NewStr);
	static void ParseEStr(const string &EStr, string &Mantissa, string &Exponent);
	static bool AllNines(const string &s);
	static bool AllZeros(const string &s);
	static bool OneZeros(const string &s);
	static double rr(double lo, double hi);
	};

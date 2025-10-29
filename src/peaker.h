#pragma once

#include "varspec.h"

typedef double (*PTR_EVAL_FUNC)(const vector<double> &xv);

class Peaker
	{
public:
	static FILE *m_fTsv;

public:
	const Peaker *m_Parent = 0;
	const string *m_ptrWhoAmI = 0;
	vector<const VarSpec *> m_VarSpecs;
	double m_Target_dy = DBL_MAX;
	double m_Min_dy = DBL_MAX;
	double m_Max_dy = DBL_MAX;
	double m_Min_Height = DBL_MAX;
	uint m_SigFig = 2;
	uint m_LatinBinCount = 0;
	uint m_EvaluateCacheHits = 0;
	uint m_RequestIdx = 0;
	double m_Change = 1.5;
	double m_Noise = 0.1;
	PTR_EVAL_FUNC m_EvalFunc = 0;
	//string m_Msg;
	string m_Cmd;
	string m_QueueDir;
	bool m_SkipInit = false;

	vector<vector<double> > m_xvs;		// Values tried
	vector<double> m_ys;				// Results
	vector<string> m_whys;				// Reason for Eval
	double m_Best_y = DBL_MAX;
	vector<double> m_Best_xv;
	time_t m_LastImprovedTime = 0;
	uint m_LastImprovedEvalCount = 0;

	vector<double> m_Deltas;

// Hooke-Jeeves parameters
	uint m_HJ_MaxExtendIters = 100;
	uint m_HJ_MaxIters = 100;

// Hooke-Jeeves state
	uint m_HJ_Direction = UINT_MAX;		// current axis
	bool m_HJ_ExtendPlus = false;	// +/- direction

public:
	Peaker(const Peaker *Parent, const string *WhoAmI)
		{
		m_Parent = Parent;
		m_ptrWhoAmI = WhoAmI;
		}

	void Clear()
		{
		m_VarSpecs.clear();
		m_Target_dy = DBL_MAX;
		m_Min_dy = DBL_MAX;
		m_Max_dy = DBL_MAX;
		m_Min_Height = DBL_MAX;
		m_SigFig = 2;
		m_SigFig = UINT_MAX;
		m_LatinBinCount = 0;
		m_EvaluateCacheHits = 0;
		m_RequestIdx = 0;
		m_EvalFunc = 0;
		m_Cmd.clear();
		m_QueueDir.clear();

		m_xvs.clear();
		m_ys.clear();
		m_whys.clear();
		m_Best_y = DBL_MAX;
		m_LastImprovedTime = 0;
		m_LastImprovedEvalCount = 0;
		m_Best_xv.clear();

		m_Deltas.clear();
		}

public:
	void Init(const vector<string> &SpecLines, PTR_EVAL_FUNC EF);
	void WriteStatusPage(FILE *f) const;
	void LogSpecs() const;
	void InitDeltas();
	void LogDeltas() const;
	void CleanQueue();
	uint GetSigFig(uint VarIdx) const;
	const char *VarsToStr(const vector<double> &xv, string &s,
						  const string sep=";") const;
	const char *VarToStr(double x, uint VarIdx, string &s) const;
	const char *GetVarName(uint VarIdx) const;
	bool Cmp_xs(const vector<double> &xs1, const vector<double> &xs2) const;
	void Round_xs(const vector<double> &xv, vector<double> &Rounded_xv) const;
	uint Find_xs(const vector<double> &xv) const;
	uint GetVarCount() const { return SIZE(m_VarSpecs); }
	const VarSpec &GetVarSpec(uint VarIdx) const;
	void GetLatinHypercube(vector<vector<double> > &xvs);
	double Evaluate(const vector<double> &xv, const string &why);
	double Calc(const vector<double> &xv);
	double CalcQueue(const vector<double> &xv);
	void Run();
	void LogState() const;
	void GetBestVarStr(string &s) const;
	void GetBestVars(vector<double> &xs) const;
	double GetRounded(double x, uint SigFig) const;
	void RunInitialValues();
	void RunLatin(vector<vector<double> > &xvs, vector<double> &ys);
	double IncreaseWithNoise() const;
	double rr(double lo, double hi) const;
	bool VarIsConstant(uint VarIdx) const;
	double GetMinDelta(uint VarIdx) const;
	double GetLatinBinWidth(uint VarIdx) const;
	double GetLatinValueByBinIdx(uint VarIdx, uint BinIdx, uint BinCount) const;
	void LogLatinBins() const;
	void RunNestedLatin(uint TopN);
	void GetPeakerPathStr(string &s) const;
	void AppendChildResult(const vector<double> &xv, double y,
		const string &childname, const string &why);
	bool AdjustDeltaUp(uint VarIdx);
	bool AdjustDeltaDown(uint VarIdx);

// Hypercross
	void HC_RunHyperCross();
	uint HC_MakeHyperCross(uint xIdx, vector<double> &DeltasUp,
	  vector<double> &DeltasDown, vector<uint> &xIdx_Us, vector<uint> &xIdx_Ds);
	void HC_SearchAxis(uint VarIdx, double &DeltaUp, double &DeltaDown,
	  uint &xIdx_D, uint &xIdx_C, uint &xIdx_U);
	double DeltaVar(uint VarIdx, bool Plus, double OldValue);

// Hooke-Jeeves
	void HJ_RunHookeJeeves();
	void HJ_Explore();
	void HJ_Extend();
	double HJ_TryDelta(const string &reason,
		const vector<double> &Start_xv, uint VarIdx, bool Plus);

public:
	static void GetLatinHypercubeIdxs(uint VarCount, uint BinCount,
		vector<vector<uint> > &IdxMx);
	};

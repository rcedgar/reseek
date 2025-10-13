#pragma once

#include "varspec.h"

typedef double (*PTR_EVAL_FUNC)(const vector<double> &xv);

class Peaker
	{
public:
	static FILE *m_fTsv;

public:
	vector<const VarSpec *> m_VarSpecs;
	double m_Target_dy = DBL_MAX;
	double m_Min_dy = DBL_MAX;
	double m_Max_dy = DBL_MAX;
	double m_Min_Height = DBL_MAX;
	uint m_DefaultSigFig = 6;
	uint m_SigFig = UINT_MAX;
	bool m_DoLatin = false;
	uint m_EvaluateCacheHits = 0;
	uint m_RequestIdx = 0;
	double m_Change = 1.5;
	double m_Noise = 0.1;
	PTR_EVAL_FUNC m_EvalFunc = 0;
	string m_Msg;
	string m_Cmd;
	string m_QueueDir;

	vector<vector<double> > m_xvs;		// Values tried
	vector<double> m_ys;				// Results
	double m_Best_y = DBL_MAX;

	vector<double> m_Deltas;

	uint m_EvalIdx = 0;
	uint m_Best_xIdx = UINT_MAX;
	uint m_LastImprovedEvalIdx = 0;
	vector<double> m_Best_ys;
	vector<uint> m_Best_EvalIdxs;
	vector<uint> m_Best_xIdxs;

	time_t m_LastProgressTime = 0;

// Hooke-Jeeves parameters
	uint m_HJ_MaxExtendIters = 100;
	uint m_HJ_MaxIters = 100;

// Hooke-Jeeves state
	uint m_HJ_Direction = UINT_MAX;	// current axis
	bool m_HJ_Plus = false;			// +/- direction

	void Clear()
		{
		m_VarSpecs.clear();
		m_Target_dy = DBL_MAX;
		m_Min_dy = DBL_MAX;
		m_Max_dy = DBL_MAX;
		m_Min_Height = DBL_MAX;
		m_DefaultSigFig = 6;
		m_SigFig = UINT_MAX;
		m_DoLatin = false;
		m_EvaluateCacheHits = 0;
		m_RequestIdx = 0;
		m_EvalFunc = 0;
		m_Cmd.clear();
		m_QueueDir.clear();

		m_xvs.clear();
		m_ys.clear();
		m_Best_y = DBL_MAX;

		m_Deltas.clear();

		m_EvalIdx = 0;
		m_Best_xIdx = UINT_MAX;
		m_LastImprovedEvalIdx = 0;
		m_Best_ys.clear();
		m_Best_EvalIdxs.clear();
		m_Best_xIdxs.clear();

		m_LastProgressTime = 0;
		m_Msg.clear();
		}

public:
	void Init(const vector<string> &SpecLines, PTR_EVAL_FUNC EF);

	void InitDeltas();
	void LogDeltas() const;
	void CleanQueue();
	uint GetSigFig(uint VarIdx) const;
	const char *VarsToStr(const vector<double> &xv, string &s,
						  const string sep=",") const;
	const char *VarToStr(double x, uint VarIdx, string &s) const;
	const char *GetVarName(uint VarIdx) const;
	bool Cmp_xs(const vector<double> &xs1, const vector<double> &xs2) const;
	uint Find_xs(const vector<double> &xv) const;
	uint Add_xv(const vector<double> &xv);
	const vector<double> &Get_xv(uint xIdx) const;
	double Get_x(uint xIdx, uint VarIdx) const;
	double Get_y(uint xIdx) const;
	uint GetVarCount() const { return SIZE(m_VarSpecs); }
	const VarSpec &GetVarSpec(uint VarIdx) const;
	void GetLatinHypercube(vector<vector<double> > &xvs);
	bool GetLatinHypercubeVar(vector<double> &xv,
	  vector<vector<bool> > &FilledMx, uint VarIdx);
	double Evaluate(uint xIdx, bool UnsetOk = false);
	double Calc(const vector<double> &xv);
	double CalcQueue(const vector<double> &xv);
	void Run();
	void LogPair(uint xIdx1, uint xIdx2) const;
	void LogState() const;
	void ProgressLogSummary() const;
	void GetBestVarStr(string &s) const;
	double GetRounded(double x, uint SigFig) const;
	void RunInitialValues();
	void RunLatin();
	double ChangeWithNoise() const;
	double rr(double lo, double hi) const;

// Hypercross
	void HC_RunHyperCross();
	uint HC_MakeHyperCross(uint xIdx, vector<double> &DeltasUp,
	  vector<double> &DeltasDown, vector<uint> &xIdx_Us, vector<uint> &xIdx_Ds);
	void HC_SearchAxis(uint VarIdx, double &DeltaUp, double &DeltaDown,
	  uint &xIdx_D, uint &xIdx_C, uint &xIdx_U);

// Hooke-Jeeves
	void HJ_RunHookeJeeves();
	double HJ_Explore();
	void HJ_Extend();
	double HJ_EvalDelta(uint xIdx, uint VarIdx, bool Plus);
	};

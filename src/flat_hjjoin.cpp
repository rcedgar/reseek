#include "myutils.h"
#include "statsig.h"
#include "parabench.h"
#include "peaker.h"

float *read_join_data(
	const string &fn,
	const lookup &look,
	vector<string> &names);

static Peaker *s_Peaker;
static float *s_join_data;
static vector<string> s_join_names;
static vector<string> s_var_names;
static FastBench *s_FB;

static uint s_varidx_mega = UINT_MAX;
static uint s_varidx_megarev = UINT_MAX;
static uint s_varidx_megaselfrev = UINT_MAX;
static uint s_varidx_entropy = UINT_MAX;
static uint s_varidx_lddt = UINT_MAX;
static uint s_varidx_dali = UINT_MAX;
static uint s_varidx_dalix = UINT_MAX;
static uint s_varidx_lpow = UINT_MAX;
static uint s_varidx_ladd = UINT_MAX;
static uint s_varidx_dpw = UINT_MAX;
static uint s_varidx_revtsw = UINT_MAX;
static uint s_varidx_lddtw = UINT_MAX;
static uint s_varidx_logladd = UINT_MAX;

static uint s_joinidx_q_selfrev_mega = UINT_MAX;
static uint s_joinidx_t_selfrev_mega = UINT_MAX;
static uint s_joinidx_q_selfrev_nu = UINT_MAX;
static uint s_joinidx_t_selfrev_nu = UINT_MAX;
static uint s_joinidx_q_L = UINT_MAX;
static uint s_joinidx_t_L = UINT_MAX;

static uint s_joinidx_mega = UINT_MAX;
static uint s_joinidx_megarev = UINT_MAX;
static uint s_joinidx_nu = UINT_MAX;
static uint s_joinidx_dali = UINT_MAX;
static uint s_joinidx_dalix = UINT_MAX;
static uint s_joinidx_lddt = UINT_MAX;
static uint s_joinidx_entropy = UINT_MAX;
static uint s_joinidx_l2 = UINT_MAX;
static uint s_joinidx_l = UINT_MAX;
static uint s_joinidx_dpscore = UINT_MAX;
static uint s_joinidx_selfrev = UINT_MAX;
static uint s_joinidx_newts = UINT_MAX;

static double calc_score(
	const float *data,
	vector<double> &values)
	{
	const uint nvar = uint(s_var_names.size());
	asserta(values.size() == nvar);

#define x(nm)	float value_##nm = (s_varidx_##nm == UINT_MAX ? 0 : (float) values[s_varidx_##nm])
	x(mega);
	x(megarev);
	x(megaselfrev);
	x(dalix);
	x(entropy);
	x(lddt);
#undef x

//[ 0]  q_selfrev_mega
//[ 1]  t_selfrev_mega
//[ 2]  q_selfrev_nu
//[ 3]  t_selfrev_nu
//[ 4]  q_L
//[ 5]  t_L
//[ 6]  mega
//[ 7]  megarev
//[ 8]  nu
//[ 9]  dali
//[10]  dalix
//[11]  entropy
//[12]  lddt
#define x(nm)	float nm = (s_joinidx_##nm == UINT_MAX ? 0 : data[s_joinidx_##nm]);
	x(q_selfrev_mega);
	x(t_selfrev_mega);
	x(q_selfrev_nu);
	x(t_selfrev_nu);
	x(dpscore);
	x(mega);
	x(megarev);
	x(nu);
	x(entropy);
	x(lddt);
	x(dali);
	x(dalix);
	x(q_L);
	x(t_L);
#undef x

	asserta(s_joinidx_q_selfrev_mega != UINT_MAX);
	asserta(s_joinidx_t_selfrev_mega != UINT_MAX);
	asserta(s_joinidx_mega != UINT_MAX);
	asserta(s_joinidx_megarev != UINT_MAX);
	asserta(s_joinidx_dalix != UINT_MAX);
	asserta(s_joinidx_entropy != UINT_MAX);
	asserta(s_joinidx_lddt != UINT_MAX);

	//asserta(s_varidx_mega != UINT_MAX);
	asserta(s_varidx_megaselfrev != UINT_MAX);
	asserta(s_varidx_dalix != UINT_MAX);
	asserta(s_varidx_entropy != UINT_MAX);

	float megaselfrev = (q_selfrev_mega + t_selfrev_mega)/2;

	float score = mega;
	score -= value_megarev*megarev;
	score -= value_megaselfrev*megaselfrev;
	score += value_entropy*entropy/250;
	score += value_dalix*dalix*10;
	score += value_lddt*lddt/4;

	asserta(!isnan(score));
	return score;
	}

static double get_value(
	const vector<string> &xv,
	uint idx,
	const string &name)
	{
	asserta(s_Peaker->m_VarNames[idx] == name);
	return StrToFloat(xv[idx]);
	}

static double EvalSum3(const vector<string> &xv)
	{
	asserta(s_FB);
	vector<double> weights;
	s_Peaker->xv2values(xv, weights);

	uint nf = uint(s_join_names.size());
	const uint npair =
		s_FB->m_look->get_pair_count_upper_triangle_with_diagonal();
	bool do_dope = (s_FB->m_dope != 0);
	for (uint k = 0; k < npair; ++k)
		{
		if (do_dope && !s_FB->in_dope(k))
			{
			s_FB->m_Scores[k] = -9999;
			continue;
			}
		float score = (float) calc_score(s_join_data + nf*k, weights);
		s_FB->m_Scores[k] = score;
		}
	Progress("Sorting...\r");
	s_FB->SetScoreOrder_Parallel();
	double Sum3 = s_FB->Bench();
	return Sum3;
	}

static void Optimize(
	const vector<string> &SpecLines,
	double &Best_y,
	vector<string> &Best_xv)
	{
	const string OptName("latinclimb");
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	uint LatinBinCount = Peaker::SpecGetInt(GlobalSpec, "latin", UINT_MAX);
	uint HJCount = Peaker::SpecGetInt(GlobalSpec, "hj", UINT_MAX);
	asserta(LatinBinCount != UINT_MAX);
	asserta(HJCount != UINT_MAX);

	Peaker &P = *new Peaker(0, OptName);
	P.Init(SpecLines, EvalSum3);
	s_Peaker = &P;

	ProgressLog("=========================================\n");
	ProgressLog("%s latin (%u)\n", OptName.c_str(), LatinBinCount);
	ProgressLog("=========================================\n");
	asserta(LatinBinCount > 0);
	P.RunLatin(LatinBinCount);

	vector<uint> TopEvalIdxs;
	P.GetTopEvalIdxs(HJCount, TopEvalIdxs);
	const uint n = SIZE(TopEvalIdxs);
	if (n == 0)
		Die("No evals %s", OptName.c_str());

	for (uint k = 0; k < n; ++k)
		{
		ProgressLog("=========================================\n");
		ProgressLog("%s HJ %u/%u\n", OptName.c_str(), k+1, n);
		ProgressLog("=========================================\n");

		uint EvalIdx = TopEvalIdxs[k];
		string ChildName;
		Ps(ChildName, "HJ%u/%u", k+1, n);
		Peaker *Child = P.MakeChild(ChildName);
		double y = P.m_ys[EvalIdx];
		const vector<string> &xv = P.m_xvs[EvalIdx];
		Child->AppendResult(xv, y, "HJstart");
		Child->HJ_RunHookeJeeves();
		P.AppendChildResults(*Child);
		delete Child;
		ProgressLog("=========================================\n");
		ProgressLog("%s HJ %u/%u converged\n", OptName.c_str(), k+1, n);
		ProgressLog("=========================================\n");
		}
	Best_y = P.m_Best_y;
	Best_xv = P.m_Best_xv;
	ProgressLog("=========================================\n");
	ProgressLog("%s completed\n", OptName.c_str());
	ProgressLog("=========================================\n");
	}

static void Climb(ParaBench &PS, const vector<string> &SpecLines)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	vector<string> Fields;
	string ParamStr;
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		{
		const string &Line = SpecLines[i];
		if (StartsWith(Line, "#init "))
			{
			ParamStr = Line.substr(6);
			break;
			}
		}
	if (ParamStr.empty())
		Die("Missing #init in spec");

	vector<string> Fields2, VarNames, Init_xv;
	Split(ParamStr, Fields, ';');
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		Split(Fields[i], Fields2, '=');
		asserta(SIZE(Fields2) == 2);
		VarNames.push_back(Fields2[0]);
		Init_xv.push_back(Fields2[1]);
		}
	const uint VarCount = SIZE(VarNames);

	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	s_Peaker = &Pfull;

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

void cmd_flat_hjjoin()
	{
	asserta(optset_lookup);
	asserta(!optset_spec);

	const string &joinfn = g_Arg1;

	s_FB = new FastBench;
	s_FB->m_scores_are_evalues = opt(scores_are_evalues);
	s_FB->ReadLookup(opt(lookup));
	if (optset_dope)
		s_FB->ReadDope(opt(dope));
	s_FB->Alloc();

	s_join_data = read_join_data(joinfn, *s_FB->m_look, s_join_names);
	for (uint i = 0; i < s_join_names.size(); ++i)
		{
		const string &name = s_join_names[i];
		if (0) ;
#define x(nm)	else if (name == #nm) s_joinidx_##nm = i
		x(q_selfrev_mega);
		x(t_selfrev_mega);
		x(q_selfrev_nu);
		x(t_selfrev_nu);
		x(mega);
		x(megarev);
		x(nu);
		x(dali);
		x(dalix);
		x(lddt);
		x(entropy);
		x(q_L);
		x(t_L);
		x(l2);
		x(l);
		x(dpscore);
		x(selfrev);
		x(newts);
#undef x
		else Die("var=%s", name.c_str());
		}

	vector<string> SpecLines;
	SpecLines.push_back("strategy=latinclimb;");
	SpecLines.push_back("latin=32;");
	SpecLines.push_back("rates=1.3,1.05;");
	SpecLines.push_back("hj=1;");
	//SpecLines.push_back("var=mega;min=0;max=1;weight=no;sigfig=3;");
	SpecLines.push_back("var=megarev;min=0;max=1;weight=no;sigfig=3;");
	SpecLines.push_back("var=megaselfrev;min=0;max=1;weight=no;sigfig=3;");
	SpecLines.push_back("var=dalix;min=0;max=1;weight=no;sigfig=3;");
	SpecLines.push_back("var=entropy;min=0;max=1;weight=no;sigfig=3;");
	SpecLines.push_back("var=lddt;min=0;max=1;weight=no;sigfig=3;");

	Peaker::GetVarNames(SpecLines, s_var_names);
	for (uint i = 0; i < s_var_names.size(); ++i)
		{
		const string &name = s_var_names[i];
		if (0) ;
#define x(nm)	else if (name == #nm) s_varidx_##nm = i
		x(mega);
		x(megarev);
		x(megaselfrev);
		x(entropy);
		x(lddt);
		x(dali);
		x(dalix);
		x(lpow);
		x(ladd);
		x(dpw);
		x(revtsw);
		x(lddtw);
		x(logladd);
#undef x
		else Die("var=%s", name.c_str());
		}

	if (optset_varstr)
		{
		Peaker &P = *new Peaker(0, "");
		P.Init(SpecLines, EvalSum3);
		s_Peaker = &P;

		vector<string> xv;
		s_Peaker->xss2xv(opt(varstr), xv);
		EvalSum3(xv);
		return;
		}

	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	double Best_y;
	vector<string> Best_xv;
	Optimize(SpecLines, Best_y, Best_xv);
	string xstr;
	s_Peaker->xv2xss(Best_xv, xstr);
	Log("@TSV@");
	Log("\t%.4g", Best_y);
	Log("\t%s", xstr.c_str());
	Log("\n");
	Progress("\n");
	Progress("\n");
	ProgressLog("FINAL [%.4g]  %s\n", Best_y, xstr.c_str());
	Progress("\n");
	Progress("\n");
	vector<double> best_values;
	s_Peaker->xv2values(Best_xv, best_values);
	for (uint i = 0; i < best_values.size(); ++i)
		ProgressLog("%s=%.4g\n",
			s_Peaker->m_VarNames[i].c_str(),
			best_values[i]);
	}

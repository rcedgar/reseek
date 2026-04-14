#include "myutils.h"
#include "statsig.h"
#include "parasearch.h"
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

/***
-rw-r--r-- 1 bob bob 744 Apr 13  2026 ../2026-04-12_ts_feature_vector/feature_stats.txt

selfrev_mega  mean  14.01,  med  14.57,  min  11.22,  max  14.57
  selfrev_nu  mean  11.56,  med  11.02,  min  4.786,  max  50.82
        mega  mean  131.8,  med    115,  min    115,  max    216
     megarev  mean  103.3,  med     96,  min     36,  max    589
          nu  mean  10.16,  med   9.03,  min  0.929,  max    264
        dali  mean  7.788,  med   7.79,  min    1.3,  max   19.1
     entropy  mean  104.3,  med     93,  min      5,  max   1550
        lddt  mean  32.02,  med     26,  min   -159,  max    841
***/

static uint s_varidx_mega = UINT_MAX;
static uint s_varidx_megarev = UINT_MAX;
static uint s_varidx_megaselfrev = UINT_MAX;
static uint s_varidx_entropy = UINT_MAX;
static uint s_varidx_lddt = UINT_MAX;
static uint s_varidx_dali = UINT_MAX;
static uint s_varidx_lpow = UINT_MAX;
static uint s_varidx_ladd = UINT_MAX;

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
static uint s_joinidx_lddt = UINT_MAX;
static uint s_joinidx_entropy = UINT_MAX;

static double calc_score(
	const float *data,
	vector<double> &weights)
	{
	const uint nvar = uint(s_var_names.size());
	asserta(weights.size() == nvar);

#define x(nm)	float w_##nm = (s_varidx_##nm == UINT_MAX ? 0 : (float) weights[s_varidx_##nm])
	x(mega);
	x(megarev);
	x(megaselfrev);
	x(entropy);
	x(lddt);
	x(dali);
	x(lpow);
	x(ladd);
#undef x

#define x(nm)	float nm = (s_joinidx_##nm == UINT_MAX ? 0 : data[s_joinidx_##nm]);
	x(q_selfrev_mega);
	x(t_selfrev_mega);
	x(q_selfrev_nu);
	x(t_selfrev_nu);
	x(mega);
	x(megarev);
	x(nu);
	x(entropy);
	x(lddt);
	x(dali);
	x(q_L);
	x(t_L);
#undef x

	if (s_varidx_lpow != UINT_MAX)
		{
		asserta(s_joinidx_q_L != UINT_MAX);
		asserta(s_joinidx_t_L != UINT_MAX);
		float L = (q_L + t_L)/2;
		float score = mega*pow(L, w_lpow);
		return score;
		}

	if (s_varidx_ladd != UINT_MAX &&
		s_varidx_mega == UINT_MAX &&
		s_varidx_megarev == UINT_MAX)
		{
		asserta(s_joinidx_q_L != UINT_MAX);
		asserta(s_joinidx_t_L != UINT_MAX);
		float L = (q_L + t_L)/2;
		float score = mega/(L + w_ladd*500);
		return score;
		}

	if (s_varidx_ladd != UINT_MAX &&
		s_varidx_mega != UINT_MAX &&
		s_varidx_megarev != UINT_MAX)
		{
		asserta(s_joinidx_q_L != UINT_MAX);
		asserta(s_joinidx_t_L != UINT_MAX);
		float L = (q_L + t_L)/2;
		float score = (w_mega*mega - w_megarev*megarev)/(L + w_ladd*500);
		return score;
		}

	float megaselfrev = (q_selfrev_mega + t_selfrev_mega)/2;
	//float nurev = (q_selfrev_nu + t_selfrev_nu)/2;

	float score = 0;
	score += w_mega*mega*20;
	//score += w_nu*nu/3;
	score += w_dali*dali/5;
	score += w_entropy*entropy/5;
	score += w_lddt*lddt*20;
	score -= w_megarev*megarev*20;
	score -= w_megaselfrev*megaselfrev*20;
	//score -= w_nurev*nurev*20;
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
	s_FB->Bench();
	float Sum3 = s_FB->m_Sum3;
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

static void Climb(ParaSearch &PS, const vector<string> &SpecLines)
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
	asserta(optset_spec);

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
		x(lddt);
		x(entropy);
		x(q_L);
		x(t_L);
#undef x
		else Die("var=%s", name.c_str());
		}

	vector<string> SpecLines;
	ReadLinesFromFile(opt(spec), SpecLines);

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
		x(lpow);
		x(ladd);
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
	Progress("FINAL [%.4g]  %s", Best_y, xstr.c_str());
	Progress("\n");
	}

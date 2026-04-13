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
static FastBench *s_FB;

/***
                mean_tp  mean_fp   std_tp  std_fp
feature                                                                                 
q_selfrev_mega  10.9941  11.5166   2.6009  3.2537	[0]
t_selfrev_mega  10.9874  11.1899   2.6918  3.1715	[1]
q_selfrev_nu   103.8626 103.1983  39.1769 37.6556	[2]
t_selfrev_nu   102.4237 103.7084  36.6551 38.6891	[3]
mega            33.4670   9.1704  35.2132  2.8570	[4]
nu             298.1936  94.5640 272.5236 30.5869	[5]
dali           151.3461  20.2054 208.7387 22.2518	[6]
entropy        284.3277  51.2113 233.6193 57.6467	[7]
lddt             7.0773   1.8915   5.3000  0.8488	[8]
***/

static double calc_score(
	const float *data,
	vector<double> &weights)
	{
	asserta(weights.size() == 5);

	//float w_nu = (float) weights[1];
	//float w_nurev = (float) weights[6];

	float w_mega = (float) weights[0];
	float w_dali = (float) weights[1];
	float w_entropy = (float) weights[2];
	float w_lddt = (float) weights[3];
	float w_megarev = (float) weights[4];

	float q_selfrev_mega = data[0];
	float t_selfrev_mega = data[1];
	//float q_selfrev_nu = data[2];
	//float t_selfrev_nu = data[3];
	float mega = data[4];
	float nu = data[5];
	float dali  = data[6];
	float entropy = data[7];
	float lddt = data[8];

	float megarev = (q_selfrev_mega + t_selfrev_mega)/2;
	//float nurev = (q_selfrev_nu + t_selfrev_nu)/2;

	float score = 0;
	score += w_mega*mega*20;
	//score += w_nu*nu/3;
	score += w_dali*dali/5;
	score += w_entropy*entropy/5;
	score += w_lddt*lddt*20;
	score -= w_megarev*megarev*20;
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
	asserta(xv.size() == 5);
	vector<double> weights(5);

	//weights[1] = get_value(xv, 1, "nu");
	//weights[6] = get_value(xv, 6, "nurev");

	weights[0] = get_value(xv, 0, "mega");
	weights[1] = get_value(xv, 1, "dali");
	weights[2] = get_value(xv, 2, "entropy");
	weights[3] = get_value(xv, 3, "lddt");
	weights[4] = get_value(xv, 4, "megarev");

	uint nf = uint(s_join_names.size());
	const uint npair =
		s_FB->m_look->get_pair_count_upper_triangle_with_diagonal();
	for (uint i = 0; i < npair; ++i)
		{
		float score = (float) calc_score(s_join_data + nf*i, weights);
		s_FB->m_Scores[i] = score;
		}
	ProgressLog("Sorting...");
	s_FB->SetScoreOrder_Parallel();
	ProgressLog("\n");
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
	asserta(optset_lookup);

	const string &joinfn = g_Arg1;

	s_FB = new FastBench;
	s_FB->m_scores_are_evalues = opt(scores_are_evalues);
	s_FB->ReadLookup(opt(lookup));
	s_FB->Alloc();

	s_join_data = read_join_data(joinfn, *s_FB->m_look, s_join_names);

	vector<string> SpecLines;
	SpecLines.push_back("strategy=latinclimb;");
	SpecLines.push_back("latin=16;");
	SpecLines.push_back("rates=1.3,1.05;");
	SpecLines.push_back("hj=1;");
	SpecLines.push_back("var=mega;min=0;max=1;weight=yes;sigfig=3;");
	SpecLines.push_back("var=dali;min=0;max=1;weight=yes;sigfig=3;");
	SpecLines.push_back("var=entropy;min=0;max=1;weight=yes;sigfig=3;");
	SpecLines.push_back("var=lddt;min=0;max=1;weight=yes;sigfig=3;");
	SpecLines.push_back("var=megarev;min=0;max=1;weight=yes;sigfig=3;");

	double Best_y;
	vector<string> Best_xv;
	Optimize(SpecLines, Best_y, Best_xv);
	}

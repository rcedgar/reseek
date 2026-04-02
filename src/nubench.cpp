#include "myutils.h"
#include "statsig.h"
#include "parasearch.h"
#include "flat_bench.h"

void GetFeatures(
	const string &varstr,
	vector<string> &feature_names,
	vector<float> &weights)
	{
	void ParseVarStr(
		const string &VarStr,
		vector<string> &Names,
		vector<float> &Values);

	vector<string> names;
	vector<float> values;
	ParseVarStr(varstr, names, values);
	vector<string> AlphaNames;

	vector<string> scalar_names;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(
		names, values,
		feature_names, weights,
		scalar_names, scalar_values);
	}

void cmd_numx()
	{
	vector<string> feature_names;
	vector<float> weights;
	GetFeatures(g_Arg1, feature_names, weights);
	const uint nfeat = uint(feature_names.size());
	asserta(nfeat > 0);
	asserta(weights.size() == nfeat);

	flat_features ff;
	ff.init(feature_names);
	asserta(ff.m_nfeat == nfeat);

	unordered_map<string, float> name2weight;
	for (uint fi = 0; fi < ff.m_nfeat; ++fi)
		name2weight[feature_names[fi]] = weights[fi];

	Paralign::set_flat_compound(ff, name2weight, 1, 1, 1, 1);
	Paralign::LogMatrix();
	Paralign::LogSWFastMatrix();
	}

void cmd_nubench()
	{
	asserta(optset_mxpattern);
	vector<string> feature_names;
	vector<float> weights;
	GetFeatures(g_Arg1, feature_names, weights);
	const uint nfeat = uint(feature_names.size());
	asserta(nfeat > 0);
	asserta(weights.size() == nfeat);
	asserta(!optset_scale);

	flat_features &ff = ParaSearch::m_ff;
	ff.init(feature_names);
	asserta(ff.m_nfeat == nfeat);
	ff.read_logoddsvec_pattern(opt(mxpattern));

	unordered_map<string, float> name2weight;
	for (uint fi = 0; fi < ff.m_nfeat; ++fi)
		name2weight[feature_names[fi]] = weights[fi];

	string AlignMethod = "para";
	if (optset_alignmethod)
		AlignMethod = string(opt(alignmethod));

	asserta(optset_db);
	const string &DBFN = opt(db);

	float Scale = 1.0f;
	int IntOpen = 2;
	int IntExt = 1;
	if (optset_scalef) Scale = float(opt(scalef));
	if (optset_intopen) IntOpen = opt(intopen);
	if (optset_intext) IntExt = opt(intext);
	Paralign::set_flat_compound(ff, name2weight, Scale, 
		IntOpen, IntExt, 777);

	ParaSearch PS;
	PS.GetByteSeqs(DBFN, "nuletters");
	PS.SetLookupFromLabels();
	PS.Search(AlignMethod, false);
	PS.SetScoreOrder();
	PS.WriteHits(opt(output));
	PS.Bench();
	}

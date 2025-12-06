#include "myutils.h"
#include "statsig.h"
#include "parasearch.h"

static void GetFeatures(const string &s,
	vector<FEATURE> &Fs, vector<float> &Weights)
	{
	Fs.clear();
	Weights.clear();
	vector<string> Fields;
	Split(s, Fields, ';');
	const uint n = SIZE(Fields);
	for (uint Idx = 0; Idx < n; ++Idx)
		{
		vector<string> Fields2;
		Split(Fields[Idx], Fields2, '=');
		asserta(SIZE(Fields2) == 2);
		FEATURE F = StrToFeature(Fields2[0].c_str());
		float Weight = StrToFloatf(Fields2[1]);
		Fs.push_back(F);
		Weights.push_back(Weight);
		}
	}

void cmd_numx()
	{
	vector<FEATURE> Fs;
	vector<float> Weights;
	GetFeatures(g_Arg1, Fs, Weights);
	const uint FeatureCount = SIZE(Fs);
	asserta(FeatureCount > 0);
	Paralign::SetCompoundMx(Fs, Weights, 1, 1, 1, 1);
	Paralign::LogMatrix();
	Paralign::LogSWFastMatrix();
	}

void cmd_nubench()
	{
	vector<FEATURE> Fs;
	vector<float> Weights;
	GetFeatures(g_Arg1, Fs, Weights);
	const uint FeatureCount = SIZE(Fs);
	asserta(FeatureCount > 0);

	string AlignMethod = "para";
	if (optset_alignmethod)
		AlignMethod = string(opt(alignmethod));

	asserta(optset_db);
	const string &DBFN = opt(db);

	int Scale = 1;
	int IntOpen = 2;
	int IntExt = 1;
	if (optset_scale) Scale = opt(scale);
	if (optset_intopen) IntOpen = opt(intopen);
	if (optset_intext) IntExt = opt(intext);
	Paralign::SetCompoundMx(Fs, Weights, Scale, 
		IntOpen, IntExt, 777);

	ParaSearch PS;
	PS.m_NuFs = Fs;
	PS.ReadLookup(opt(lookup));
	PS.GetByteSeqs(DBFN, "nuletters");
	PS.SetDomIdxs();
	PS.Search(AlignMethod);
	PS.WriteHits(opt(output));
	PS.Bench();
	}

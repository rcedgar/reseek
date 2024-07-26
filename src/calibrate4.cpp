#include "myutils.h"
#include "scop40bench.h"
#include <unordered_map>

/***
-calibrate4
Input:
	-calibrate4 /z/int/reseek_calibrate/scop40_with_ts/hits_with_ts.tsv
	-input2     /z/int/reseek_calibrate/features_all_scop40/output.tsv
	-model      /z/a/res/reseekml/torchout/reseekml_v5.model.tsv

Output:
	/z/int/reseek_calibrate/scop40_with_ts/hits_with_model_prob.tsv
***/

void ReadNN(const string &FN);
void GetFeatures(const vector<uint> &Bins, vector<float> &Features);
double GetNNP(const vector<double> &x);

void cmd_calibrate4()
	{
	DSSParams Params;
	Params.SetFromCmdLine(10000);

	asserta(optset_model);
	asserta(optset_input2);
	asserta(optset_output);

	const float m = 20.5f;
	const float b = 2.9f;

	ReadNN(opt_model);

/***
              TS  0.0312  0.0938  0.156  0.219  0.281  0.344  0.406  0.469  0.531  0.594  0.656  0.719  0.781  0.844  0.906  0.969
d3nfka_/b.36.1.1     247     598     22      2      2      3      4      8     14     10     16      1      1      0      0      0
d1t6ca2/c.55.1.8     887     421     35     13      1      0      0      0      0      0      0      0      1      0      0      0
d2gtlm1/b.61.7.1     167     419    196     38     10      2      2      0      0      0      0      0      1      0      0      1
***/
	unordered_map<string, vector<float> > Label2Features;
	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(opt_input2);
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	asserta(StartsWith(Line, "TS\t"));
	while (ReadLineStdioFile(f, Line))
		{
		vector<uint> Bins;
		vector<float> Features;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 17);
		const string &Label = Fields[0];
		for (uint i = 0; i < 16; ++i)
			Bins.push_back(StrToUint(Fields[i+1]));
		GetFeatures(Bins, Features);
		asserta(SIZE(Features) == 16);
		Label2Features[Label] = Features;
		}
	CloseStdioFile(f);
	f = OpenStdioFile(g_Arg1);
	FILE *fOut = CreateStdioFile(opt_output);
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 4);
		const string &Query = Fields[0];
		const string &Target = Fields[1];
		//float Evalue = StrToFloatf(Fields[2]);
		double TS = StrToFloat(Fields[3]);
		unordered_map<string, vector<float> >::const_iterator iter =
		  Label2Features.find(Target);
		if (iter == Label2Features.end())
			Die("Label not found >%s", Target.c_str());
		const vector<float> &Features = iter->second;
		asserta(SIZE(Features) == 16);
		vector<double> x;
		x.push_back(TS);
		for (uint i = 0; i < 16; ++i)
			x.push_back(double(Features[i]));
		double P = GetNNP(x);
		double E = 10000*(1 - P)*1e-6;
		double Ed = Params.GetEvalue((float) TS);
		if (Ed > 1.5) //works quite well
			E = Ed;
		fprintf(fOut, "%s\t%s\t%.3g\t%.3g\n", Query.c_str(), Target.c_str(), E, Ed);
		}
	CloseStdioFile(fOut);
	}

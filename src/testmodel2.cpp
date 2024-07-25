#include "myutils.h"

static vector<vector<double> > linear_one_weight(2);
static vector<double> linear_one_bias(2);
static uint feature_count;

// y = x W^T + b
//  W^T = transpose of weight matrix
static void forward(const vector<double> &x, vector<double> &yhat)
	{
	yhat.clear();
	asserta(SIZE(x) == feature_count);
	for (uint i = 0; i < 2; ++i)
		{
		double Sum = 0;
		for (uint j = 0; j < feature_count; ++j)
			Sum += x[j]*linear_one_weight[i][j];
		Sum += linear_one_bias[i];
		yhat.push_back(Sum);
		}
	}

static void softmax(const vector<double> &x, vector<double> &y)
	{
	y.clear();
	const uint n = SIZE(x);
	double Sum = 0;
	for (uint i = 0; i < n; ++i)
		Sum += exp(x[i]);
	if (Sum == 0)
		Sum = 1;
	for (uint i = 0; i < n; ++i)
		y.push_back(exp(x[i])/Sum);
	}

static void log_softmax(const vector<double> &x, vector<double> &y)
	{
	y.clear();
	const uint n = SIZE(x);
	double Sum = 0;
	for (uint i = 0; i < n; ++i)
		Sum += exp(x[i]);
	if (Sum == 0)
		Sum = 1;
	for (uint i = 0; i < n; ++i)
		y.push_back(x[i] - log(Sum));
	}

/***
linear_one_weight       2       17
0       9.0769  -1.0425 0.4151  -1.4175 -2.2356 -2.2793 0.2514  0.9436  0.8387  1.4274  2.7783  -2.0859 1.3354  1.0119     9.6130  -10.6313        -38.6946
1       -8.9529 1.1474  -0.3866 1.4410  1.8608  2.3841  -0.4059 -0.8979 -0.7933 -1.3125 -2.5792 1.8608  -1.0964 -1.0386    -9.4006 10.7695 38.5087
linear_one_bias 2       -2.9990 2.6514
***/

static void ReadModel(const string &FN)
	{
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 3);
	asserta(Fields[0] == "linear_one_weight");
	asserta(Fields[1] == "2");
	feature_count = StrToUint(Fields[2]);

	for (uint i = 0; i < 2; ++i)
		{
		Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == feature_count+1);
		asserta(StrToUint(Fields[0]) == i);
		linear_one_weight[i].resize(feature_count);
		for (uint j = 0; j < feature_count; ++j)
			linear_one_weight[i][j] = StrToFloat(Fields[j+1]);
		}
	Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 4);
	asserta(Fields[0] == "linear_one_bias");
	asserta(Fields[1] == "2");
	linear_one_bias[0] = StrToFloat(Fields[2]);
	linear_one_bias[1] = StrToFloat(Fields[3]);

	CloseStdioFile(f);
	}

void cmd_testmodel2()
	{
	ReadModel(g_Arg1);
	FILE *f = OpenStdioFile(opt_input);
	string Line;
	vector<string> Fields;
	for (;;)
		{
/***
			Label  T      TS     F0      F1      F2       F3       F4       F5  F6       F7  F8  F9  F10 >
d1w36b2/c.37.1.19  1  0.0224  0.998  0.0098       0  0.00245        0  0.00245   0        0   0   0    0 >
 d1kk8a2/c.37.1.9  1  0.0311  0.998       0       0        0        0        0   0        0   0   0    0 >
 d1eu8a_/c.94.1.1  1  0.0416  0.997   0.079  0.0136  0.00272  0.00545  0.00272   0  0.00272   0   0    0 >
***/
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			break;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == feature_count + 2);
		const string &Label = Fields[0];
		const string &T = Fields[1];
		if (Label == "Label")
			{
			asserta(T == "T" && Fields[2] == "TS");
			continue;
			}

		vector<double> x;
		for (uint j = 0; j < feature_count; ++j)
			x.push_back(StrToFloat(Fields[j+2]));

		vector<double> yhat;
		forward(x, yhat);

		vector<double> ysm;
		softmax(yhat, ysm);

		//vector<double> ylsm;
		//log_softmax(yhat, ylsm);

		Log("%18.18s", Label.c_str());
		Log("  y=%s", T.c_str());
		Log("  yhat=(%.3g, %.3g)", yhat[0], yhat[1]);
		Log("  ysm=(%.3g, %.3g)", ysm[0], ysm[1]);
		//Log("  ylsm=(%.3g, %.3g)", ylsm[0], ylsm[1]);
		Log("  x =");
		for (uint j = 0; j < feature_count; ++j)
			Log(" %.3g", x[j]);
		Log("\n");
		}
	CloseStdioFile(f);
	}

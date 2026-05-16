#include "myutils.h"
#include "flat_helpers.h"
#include "flat_alphas.h"
#include "flat_params.h"

/***
[scalar] 0.000435  dali
[scalar]    0.894  gap2
[scalar]   0.0466  lddt
[scalar]    0.672  revw
[scalar]    0.962  selfw

   ppack32    0.002  o
    pmdd32    0.008  o
 mendist32    0.009  o
 pendist32    0.014  ■
       pm2    0.014  ■
  nensec32    0.016  ■
 rendist32    0.022  ■
   turnd32    0.035  ■■
     sec32    0.060  ■■■■
   angle32    0.087  ■■■■■■
    pack32    0.089  ■■■■■■■
 nendist32    0.116  ■■■■■■■■■
      aa20    0.529  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
     Total    1.000

float flat_params::m_nu_filter_self_w = 0.5f;
float flat_params::m_nu_filter_rev_w = 0.27f;
int flat_params::m_nu_filter_min_fwd_score = 120;
int flat_params::m_nu_filter_min_combined_score = 43;
***/

static const string default_alpha_weights =
	"aa20=5.29E-01;"
	"angle32=8.66E-02;"
	"mendist32=9.01E-03;"
	"nendist32=1.16E-01;"
	"nensec32=1.58E-02;"
	"pack32=8.87E-02;"
	"pendist32=1.37E-02;"
	"pm2=1.40E-02;"
	"pmdd32=8.04E-03;"
	"ppack32=1.77E-03;"
	"rendist32=2.25E-02;"
	"sec32=5.97E-02;"
	"turnd32=3.49E-02;";

static const string default_gaps =
	"gap2=8.94E-01;";

static const string default_test_statistic_weights =
	"dali=4.35E-04;"
	"lddt=4.66E-02;"
	"revw=6.72E-01;"
	"selfw=9.62E-01;";

static const string default_mega_filter =
	"minfwd=15;";

static const string default_nu_filter =
	"nfselfw=0.5;"
	"nfrevw=0.27;"
	"nfminfwd=120;"
	"nfmincmb=40;";

static const string default_varstr =
	default_alpha_weights + 
	default_gaps +
	default_test_statistic_weights +
	default_mega_filter +
	default_nu_filter;

// C:\src\reseek_tune2\bash\flat_bench2_mufilter_sweep.bash
// Sum3   Secs  
// 1.761    10  minmufwd120.minmucmb60.minmgfwd20
// 1.792	23	minmufwd120.minmucmb43.minmgfwd0			
// 1.796	25	minmufwd120.minmucmb40.minmgfwd15			

void parse_varstr(
	const string &arg_VarStr,
	vector<string> &Names,
	vector<float> &Values)
	{
	Names.clear();
	Values.clear();

	string VarStr;
	if (arg_VarStr == "")
		VarStr = default_varstr;
	else
		{
		if (StartsWith(arg_VarStr, "@"))
			{
			const string fn = arg_VarStr.substr(1);
			vector<string> lines;
			ReadLinesFromFile(fn, lines);
			for (auto line : lines)
				{
				if (StartsWith(line, "#"))
					continue;
				StripWhiteSpace(line);
				asserta(EndsWith(line, ";"));
				VarStr += line;
				}
			}
		else
			VarStr = arg_VarStr;
		}

	StripAllWhiteSpace(VarStr);

	vector<string> Fields;
	Split(VarStr, Fields, ';');

	const uint n = SIZE(Fields);
	for (uint i = 0; i < n; ++i)
		{
		const string &NameEqValue = Fields[i];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			Die("SubsetBench::parse_varstr(%s) not name=value '%s'",
				VarStr.c_str(), Fields[i].c_str());
		const string &Name = Fields2[0];
		const string &ValueStr = Fields2[1];
		float Value = StrToFloatf(ValueStr);
		Names.push_back(Name);
		Values.push_back(Value);
		}
	}

void flat_classify_params(
	const vector<string> &names,
	const vector<float> &values,
	vector<string> &alpha_names,
	vector<float> &alpaha_weights,
	vector<string> &scalar_names,
	vector<float> &scalar_values)
	{
	for (uint i = 0; i < SIZE(names); ++i)
		{
		const string &name = names[i];
		float Value = values[i];
		bool is_scalar = false;

		if (name == "gap2") {is_scalar = true; }
#define x(param_name, m_name)	else if (name == #param_name) {is_scalar = true; }
#include "tunable_flat_params.h"
		if (is_scalar)
			{
			scalar_names.push_back(name);
			scalar_values.push_back(Value);
			}
		else
			{
			uint alpha_size = 0;
			FAN fan = parse_alpha_name(name, alpha_size);
			alpha_names.push_back(name);
			alpaha_weights.push_back(Value);
			}
		}
	}

void flat_make_varstr(string &varstr)
	{
	varstr.clear();

	if (feq(flat_params::m_open, flat_params::m_ext*10))
		Psa(varstr, "gap2=%.4g;\n", flat_params::m_open);
	else
		{
		Psa(varstr, "open=%.4g;\n", flat_params::m_open);
		Psa(varstr, "ext=%.4g;\n", flat_params::m_ext);
		}

#define x(param_name, member_name)	\
	if (string(#param_name) != "open" && string(#param_name) != "ext") \
		Psa(varstr, "%s=%.4g;\n", #param_name, flat_params::member_name);
#include "tunable_flat_params.h"

	for (uint fi = 0; fi < flat_alphas::m_nfeat; ++fi)
		{
		Psa(varstr, "%s=%.4g;\n",
			flat_alphas::m_alpha_names[fi],
			flat_alphas::m_weights[fi]);
		}
	}

void flat_make_peaker_spec(vector<string> &lines)
	{
	lines.clear();

	string line;
	if (feq(flat_params::m_open, flat_params::m_ext*10))
		{
		Ps(line, "var=gap2;constant=%.4g;", flat_params::m_open);
		lines.push_back(line);
		}
	else
		{
		Ps(line, "var=open;constant=%.4g;", flat_params::m_open);
		lines.push_back(line);

		Ps(line, "var=ext;constant=%.4g", flat_params::m_ext);
		lines.push_back(line);
		}

#define x(param_name, member_name)	\
	if (string(#param_name) != "open" && string(#param_name) != "ext") { \
		Ps(line, "var=%s;constant=%.4g;", #param_name, flat_params::member_name); \
		lines.push_back(line); }
#include "tunable_flat_params.h"

	for (uint fi = 0; fi < flat_alphas::m_nfeat; ++fi)
		{
		Ps(line, "var=%s;constant=%.4g;isalpha=yes;",
			flat_alphas::m_alpha_names[fi],
			flat_alphas::m_weights[fi]);
		lines.push_back(line);
		}
	}

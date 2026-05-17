#include "myutils.h"
#include "flat_helpers.h"
#include "flat_alphas.h"
#include "flat_params.h"

/***
[scalar]     0.89  gap2
[scalar]     0.95  selfw
[scalar]     0.64  revw
[scalar]    0.058  nurevw
[scalar] 0.0001755  dali
[scalar]    0.045  lddt
[scalar]    16.19  minfwd
[scalar]      0.5  nfselfw
[scalar]   0.2899  nfrevw
[scalar]    113.8  nfminfwd
[scalar]    33.82  nfmincmb

   ppack32    0.002  o
    pmdd32    0.008  o
 mendist32    0.008  o
       pm2    0.014  ■
 pendist32    0.015  ■
  nensec32    0.016  ■
 rendist32    0.022  ■
   turnd32    0.035  ■■
     sec32    0.059  ■■■■
   angle32    0.087  ■■■■■■
    pack32    0.088  ■■■■■■■
 nendist32    0.118  ■■■■■■■■■
      aa20    0.528  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
     Total    1.000
***/

static const string default_varstr =
"gap2=8.9E-01;selfw=9.5E-01;revw=6.4E-01;nurevw=5.8E-02;dali=0.00017545;lddt=4.5E-02;minfwd=16.191;nfselfw=5.0E-01;nfrevw=0.28992;nfminfwd=113.85;nfmincmb=33.82;aa20=5.2848E-01;angle32=8.6751E-02;mendist32=7.9881E-03;nendist32=1.1838E-01;nensec32=1.5953E-02;pack32=8.8342E-02;pendist32=1.4577E-02;pm2=1.3725E-02;pmdd32=7.9770E-03;ppack32=1.6439E-03;rendist32=2.1936E-02;sec32=5.9344E-02;turnd32=3.4900E-02;";

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

#include "myutils.h"
#include "flat_bench.h"
#include "flat_helpers.h"

/***
$src/reseek_tune2/bash/reduce_aa4.bash
$src/reseek_tune2/bash/aa4_final.bash

hjnumega_ACHPST-DEKNQR-FILMVWY-G.seed1.log:FINAL climb [1.23224]
intopen=2.90E+01;intext=3.00E+00;scale=8.81E+00;
aa4=5.15E-01;pm2=2.84E-01;sec32=2.00E-01;
***/

// need aa20 to construct aa4 in profiles
static const string VarStr_aa20 = "aa4=5.15E-01;pm2=2.84E-01;sec32=2.00E-01;aa20=0;";

//static const string VarStr_noaa20 = "aa4=5.15E-01;pm2=2.84E-01;sec32=2.00E-01;";

void codeseq2hexfasta(
	FILE *f, const string &label, const uint8_t *codeseq, uint L)
	{
	if (f == 0) return;
	string hexseq;
	for (uint i = 0; i < L; ++i)
		Psa(hexseq, "%02x", codeseq[i]);
	SeqToFasta(f, label, hexseq);
	}

void cmd_test_nu_codeseqs()
	{
	asserta(optset_alphadir);
	asserta(optset_lookup);
	asserta(optset_output);

	asserta(!optset_input);
	asserta(!optset_fapattern);
	asserta(!optset_mxpattern);
	asserta(!optset_spec);
	asserta(!optset_varstr);

	opt_nufilter = true;
	optset_nufilter = true;

	const string &chainsfn = g_Arg1;

	flat_bench FB;
	FB.m_nu_filter = true;
	FB.ReadLookup(opt(lookup));

	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(VarStr_aa20, param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_params params;
	params.init_from_alphadir(alphadir, alpha_names);

	vector<flat_chain_t *> chains;
	read_flat_chains(chainsfn, chains);
	FB.set_distmxs(chains);
	FB.load_profiles_chains(chains);
	FB.init_nu_filter();

	FILE *f = CreateStdioFile(opt(output));
	const uint ndom = uint(FB.m_fp.m_nu_codeseqs.size());
	asserta(ndom == FB.m_look->get_ndom());
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		const string &label = FB.m_look->get_dom(domidx);
		uint L = FB.m_fp.get_length(domidx);
		const uint8_t *codeseq = FB.m_fp.m_nu_codeseqs[domidx];
		codeseq2hexfasta(f, label, codeseq, L);
		}
	CloseStdioFile(f);
	}

void cmd_make_nu_parasail_matrix()
	{
	const string &alphadir = g_Arg1;

// intopen=2.90E+01;intext=3.00E+00;scale=8.81E+00;
// "aa4=5.15E-01;pm2=2.84E-01;sec32=2.00E-01;";
	unordered_map<string, float> name2weight;
	name2weight["aa4"] = 5.15E-01f;
	name2weight["pm2"] = 2.84E-01f;
	name2weight["sec32"] = 2.00E-01f;
	const float scale = 8.81E+00f;
	const int intopen = 29;
	const int intext = 3;
	const int saturated_score = 999;

	const vector<string> alpha_names = { "aa4", "pm2", "sec32" };

	flat_params params;
	params.init_from_alphadir(alphadir, alpha_names);

	Paralign::set_flat_compound(params, name2weight,
		scale, intopen, intext, saturated_score);
	Paralign::LogMatrix();
	}

#include "myutils.h"
#include "flat_bench.h"

/***
Component weights for compound aa4+pm2+sec32, fwd only optimized by -hjnumega
	src/2025-10_reseek_tune/2026-04-01_hjnumega_parasail/hjnumega.log
	hjnumega.log:FINAL climb [1.24933] intopen=2.30E+01;intext=3.00E+00;scale=8.39E+00;aa4=4.81E-01;pm2=3.01E-01;sec32=2.19E-01;
	=> aa4=0.481;pm2=0.301;sec32=0.219;intopen=23;intext=3;scale=8.39;
***/

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

	// aa4=0.481;pm2=0.301;sec32=0.219;intopen=23;intext=3;scale=8.39;
	const string VarStr = "aa4=0.481;pm2=0.301;sec32=0.219;";

	flat_bench FB;
	FB.m_nu_filter = true;
	FB.ReadLookup(opt(lookup));

	vector<string> param_names;
	vector<float> param_values;
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_alphas::init_from_alphadir(alphadir, alpha_names);

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

	const string VarStr = "aa4=0.481;pm2=0.301;sec32=0.219;";
	vector<string> param_names;
	vector<float> param_values;
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	flat_alphas::init_from_alphadir(alphadir, alpha_names);
	}

#include "myutils.h"
#include "chaq.h"
#include "chain_data.h"
#include "flat_params.h"
#include "flat_helpers.h"

// C:\src\py\map_nu_to_kappa.py
uint8_t g_nucode_to_kappacode[256] = {
  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,
  8,  9, 10, 11, 12, 13, 14, 15,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
  8,  9, 10, 11, 12, 13, 14, 15,  8,  9, 10, 11, 12, 13, 14, 15,
  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
 16, 17, 18, 19, 20, 21, 22, 23,  8,  9, 10, 11, 12, 13, 14, 15,
 24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
  8,  9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31,
 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23,
 24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
 24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
 24, 25, 26, 27, 28, 29, 30, 31,  8,  9, 10, 11, 12, 13, 14, 15,
 24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
};

uint8_t s_sec32_to_sec4[32] =
	{ 0, 1, 2, 3, 3, 0, 3, 0, 1, 3, 0, 3, 0, 3, 0, 3, 2, 0, 1, 0, 0, 0, 1, 1, 0, 3, 3, 0, 3, 0, 3, 1 };

static uint8_t components_to_nu(
	uint8_t aa4, uint8_t pm2, uint8_t sec32)
	{
	asserta(aa4 < 4);
	asserta(pm2 < 2);
	asserta(sec32 < 32);
	return aa4 + 4*pm2 + 4*2*sec32;
	}

static uint8_t components_to_kappa(
	uint8_t aa4, uint8_t pm2, uint8_t sec4)
	{
	asserta(aa4 < 4);
	asserta(pm2 < 2);
	asserta(sec4 < 4);
	uint kappa = aa4 + 4*pm2 + 4*2*sec4;
	asserta(kappa < 32);
	}

static void nu_to_components(
	uint8_t nu, uint8_t &aa4, uint8_t &pm2, uint8_t &sec32)
	{
	asserta(nu < 256);
	sec32 = nu/(4*2);
	pm2 = (nu - sec32*4*2)/4;
	aa4 = nu%4;
	asserta(aa4 < 4);
	asserta(pm2 < 2);
	asserta(sec32 < 32);
#if DEBUG
	uint nu_test = components_to_nu(aa4, pm2, sec32);
	assert(nu_test == nu);
#endif
	}

static uint8_t kappa_to_components(
	uint8_t kappa, uint8_t &aa4, uint8_t &pm2, uint8_t &sec4)
	{
	asserta(kappa < 32);
	sec4 = kappa/(4*2);
	pm2 = (kappa - sec4*4*2)/4;
	aa4 = kappa%4;
	asserta(aa4 < 4);
	asserta(pm2 < 2);
	asserta(sec4 < 32);
#if DEBUG
	uint kappa_test = components_to_kappa(aa4, pm2, sec4);
	assert(kappa_test == kappa);
#endif
	}

static uint8_t nu_to_kappa(uint8_t nu)
	{
	assert(nu < 256);
	uint8_t aa4, pm2, sec32;
	nu_to_components(nu, aa4, pm2, sec32);
	uint8_t sec4 = s_sec32_to_sec4[sec32];
	uint8_t kappa = aa4 + 4*pm2 + 4*2*sec4;
#if DEBUG
#endif
	return kappa;
	}

void chaq::sec32_codeseq_to_sec4(
	const uint8_t *codeseq_sec32, uint L,
	uint8_t *codeseq_sec4)
	{
	for (uint i = 0; i < L; ++i)
		{
		uint8_t sec32code = codeseq_sec32[i];
		assert(sec32code < 32);
		codeseq_sec4[i] = s_sec32_to_sec4[sec32code];
		}
	}

void chaq::codeseq_nu_to_kappa(
	const uint8_t *codeseq_nu, uint L,
	uint8_t *codeseq_kappa, size_t codeseq_kappa_bytes)
	{
	asserta(L <= codeseq_kappa_bytes);
	for (uint i = 0; i < L; ++i)
		{
		uint8_t nu_code = codeseq_nu[i];
		uint8_t kappa_code = g_nucode_to_kappacode[nu_code];
		assert(kappa_code < 32);
		codeseq_kappa[i] = kappa_code;
		}
	}

void chaq::fill_codeseq_nu(
	const char *charseq_aa20,
	const uint8_t *codeseq_pm2,
	const uint8_t *codeseq_sec32,
	const uint L,
	uint8_t *codeseq_nu,
	size_t codeseq_nu_bytes)
	{
	for (uint32_t pos = 0; pos < L; ++pos)
		{
		char c = charseq_aa20[pos];
		uint8_t code_aa20 = g_CharToLetterAmino[c];
		if (code_aa20 >= 20) code_aa20 = 0;

		const uint8_t code_pm2 = codeseq_pm2[pos];
		const uint8_t code_sec32 = codeseq_sec32[pos];

		assert(code_aa20 < 20);
		assert(code_pm2 < 2);
		assert(code_sec32 < 32);

		const uint8_t code_aa4 = chaq::m_aacode2aa4code[code_aa20];
		const uint8_t code_nu = uint8_t(code_aa4 + 4*code_pm2 + 4*2*code_sec32);

		codeseq_nu[pos] = code_nu;
		}
	}

void codeseq_to_hexfasta(FILE *f, const string &label,
	const uint8_t *codeseq, uint L)
	{
	if (f == 0) return;
	string hexseq;
	hexseq.reserve(2*L);
	for (uint pos = 0; pos < L; ++pos)
		Psa(hexseq, "%02x", codeseq[pos]);
	SeqToFasta(f, label, hexseq);
	}

void codeseq_to_fasta(FILE *f, const string &label,
	const uint8_t *codeseq, uint L, uint alpha_size)
	{
	if (f == 0) return;
	string hexseq;
	hexseq.reserve(L);
	asserta(alpha_size <= 36);
	const unsigned char *code_to_char =
		(alpha_size == 20 ? g_LetterToCharAmino : g_LetterToCharMu);

	for (uint pos = 0; pos < L; ++pos)
		hexseq += code_to_char[codeseq[pos]];
	SeqToFasta(f, label, hexseq);
	}

void cmd_test_kappa()
	{
	const string &chainfn = g_Arg1;

	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(opt(varstr), param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_classify_params(
		param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_params params;
	params.set_scalars(scalar_names, scalar_values);
	params.init_from_alphadir(alphadir, alpha_names);

	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	const uint nchain = uint(chains.size());

	FILE *f = g_fLog;
	if (f != 0)
		{
		fprintf(f, "static const uint8_t nucode_to_kappacode[256] = {\n");
		for (uint code = 0; code < 256; ++code)
			{
			uint8_t nu = uint8_t(code);
			if (nu > 0 && nu%16 == 0)
				fprintf(f, "\n");
			uint8_t kappa = nu_to_kappa(nu);
			g_nucode_to_kappacode[nu] = kappa;
			fprintf(f, " %2u,", kappa);
			}
		fprintf(f, "};");
		CloseStdioFile(f);
		}

	const uint maxL = 4000;
	const uint M = flat_params::m_distmx_bandwidth;
	sid_t *distmx = myalloc(sid_t, maxL*M);
	uint8_t *codeseq_nu = myalloc(uint8_t, maxL);
	uint8_t *codeseq_kappa = myalloc(uint8_t, maxL);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, maxL);
	FILE *fnu = CreateStdioFile("nu.hexfa");
	FILE *fkappa = CreateStdioFile("kappa.fa");
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		ProgressStep(chainidx, nchain, "Writing nu.hexfa and kappa.fa");
		const flat_chain_t *chain = chains[chainidx];
		const char *charseq_aa20 = chain->m_aa->m_data;
		const uint L = chain->get_length();
		asserta(L <= maxL);
		chaq::fill_distmx(chain, distmx);
		chaq::fill_chaq_vecs2(distmx, L, cv);
		chaq::fill_codeseq_nu(
			charseq_aa20, cv.pm2_codeseq, cv.sec32_codeseq,
			L, codeseq_nu, maxL);
		chaq::codeseq_nu_to_kappa(
			codeseq_nu, L,
			codeseq_kappa, maxL);
		codeseq_to_fasta(fkappa, chain->m_label, codeseq_kappa, L, 32);
		codeseq_to_hexfasta(fnu, chain->m_label, codeseq_nu, L);
		}
	CloseStdioFile(fnu);
	CloseStdioFile(fkappa);
	}

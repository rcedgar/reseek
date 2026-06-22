#include "myutils.h"
#include "chaq.h"
#include "chain_data.h"
#include "flat_params.h"
#include "flat_helpers.h"

uint8_t s_sec32_to_sec4[32] = {
 0, 0, 1, 1, 2, 2, 1, 3, 1, 2, 2, 2, 2, 1, 2, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 1, 2, 3, 2, 2, 2, 2
};

uint8_t g_nucode_to_kappacode[256] = {
  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,
  8,  9, 10, 11, 12, 13, 14, 15,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
  8,  9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31,
  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
 16, 17, 18, 19, 20, 21, 22, 23,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
 24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
 24, 25, 26, 27, 28, 29, 30, 31, 24, 25, 26, 27, 28, 29, 30, 31,
 24, 25, 26, 27, 28, 29, 30, 31,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
 16, 17, 18, 19, 20, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23,
};

void set_sec4_groups(const string &sec4_groups)
	{
	vector<string> flds;
	Split(sec4_groups, flds, '-');
	asserta(flds.size() == 4);
	for (uint i = 0; i < 32; ++i)
		s_sec32_to_sec4[i] = 0xff;
	for (uint code_sec4 = 0; code_sec4 < 4; ++code_sec4)
		{
		const string &fld = flds[code_sec4];
		asserta(fld.size() > 0);
		for (size_t j = 0; j < fld.size(); ++j)
			{
			char c = fld[j];
			uint8_t code_sec32 = g_CharToLetterMu[c];
			asserta(code_sec32 < 32);
			s_sec32_to_sec4[code_sec32] = code_sec4;
			}
		}
	for (uint code_sec32 = 0; code_sec32 < 32; ++code_sec32)
		asserta(s_sec32_to_sec4[code_sec32] < 32);
	}

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
	return kappa;
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

void chaq::codeseq_nu_to_kappa_inplace(
	uint8_t *codeseq, uint L)
	{
	for (uint i = 0; i < L; ++i)
		{
		uint8_t nu_code = codeseq[i];
		uint8_t kappa_code = g_nucode_to_kappacode[nu_code];
		assert(kappa_code < 32);
		codeseq[i] = kappa_code;
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

void chaq::fill_codeseq_nu_from_chain(
	const flat_chain_t *chain,
	sid_t *distmx_buffer,
	chaq_vecs2 *cv_buffer,
	uint8_t *codeseq_nu,
	uint buffer_L)
	{
	const uint L = chain->get_length();
	asserta(L <= buffer_L); // TODO=maxL

	sid_t *distmx = distmx_buffer;
	chaq_vecs2 *cv = cv_buffer;
	const char *charseq_aa20 = chain->m_aa->m_data;

	chaq::fill_distmx(chain, distmx);
	chaq::fill_chaq_vecs2(distmx, L, *cv);
	chaq::fill_codeseq_nu(
		charseq_aa20, cv->pm2_codeseq, cv->sec32_codeseq,
		L, codeseq_nu, buffer_L);
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

void cmd_kappa_fasta()
	{
	const string &chainfn = g_Arg1;
	asserta(!optset_fasta);
	if (optset_sec4_groups)
		set_sec4_groups(opt(sec4_groups));

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

	Log("uint8_t s_sec32_to_sec4[32] = { \n");
	for (uint sec32code = 0; sec32code < 32; ++sec32code)
		{
		Log(" %u", s_sec32_to_sec4[sec32code]);
		if (sec32code != 0 && sec32code != 31) Log(",");
		}
	Log("\n};\n");

	Log("uint8_t g_nucode_to_kappacode[256] = {\n");
	for (uint code = 0; code < 256; ++code)
		{
		uint8_t nu = uint8_t(code);
		if (nu > 0 && nu%16 == 0)
			Log("\n");
		uint8_t kappa = nu_to_kappa(nu);
		g_nucode_to_kappacode[nu] = kappa;
		Log(" %2u,", kappa);
		}
	Log("\n};");

	const uint M = flat_params::m_distmx_bandwidth;
	sid_t *distmx = myalloc(sid_t, flat_params::m_maxL*M);
	uint8_t *codeseq_nu = myalloc(uint8_t, flat_params::m_maxL);
	uint8_t *codeseq_kappa = myalloc(uint8_t, flat_params::m_maxL);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
	FILE *fnu = 0;
	if (optset_hexfasta)
		fnu = CreateStdioFile(opt(hexfasta));
	FILE *fkappa = 0;
	if (optset_output)
		fkappa = CreateStdioFile(opt(output));
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		ProgressStep(chainidx, nchain, "Writing nu hex / kappa fasta");
		const flat_chain_t *chain = chains[chainidx];
		const char *charseq_aa20 = chain->m_aa->m_data;
		const uint L = chain->get_length();
		asserta(L <= flat_params::m_maxL); // TODO=maxL
		chaq::fill_distmx(chain, distmx);
		chaq::fill_chaq_vecs2(distmx, L, cv);
		chaq::fill_codeseq_nu(
			charseq_aa20, cv.pm2_codeseq, cv.sec32_codeseq,
			L, codeseq_nu, flat_params::m_maxL);
		chaq::codeseq_nu_to_kappa(
			codeseq_nu, L,
			codeseq_kappa, flat_params::m_maxL);
		codeseq_to_fasta(fkappa, chain->m_label, codeseq_kappa, L, 32);
		codeseq_to_hexfasta(fnu, chain->m_label, codeseq_nu, L);
		}
	CloseStdioFile(fnu);
	CloseStdioFile(fkappa);
	}

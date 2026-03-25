#include "myutils.h"
#include "chaq.h"
#include "flat_chain.h"
#include "alpha.h"
#include "quantize.h"

static const uint M = 48;
static const uint m = 12;

static void get_values(
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	uint16_t *values)
	{
	const uint L = chain->get_length();

	uint16_t *distmx = myalloc(sid_t, L*M);
	uint16_t *pens = myalloc(uint16_t, L);
	uint16_t *mens = myalloc(uint16_t, L);
	uint16_t *nens = myalloc(uint16_t, L);
	uint16_t *rens = myalloc(uint16_t, L);
	sid_t *pensids = myalloc(sid_t, L);
	sid_t *mensids = myalloc(sid_t, L);
	sid_t *nensids = myalloc(sid_t, L);
	sid_t *rensids = myalloc(sid_t, L);

	chaq::fill_distmx(chain->m_xyz->m_data, L, M, distmx);

	chaq::fill_pen_men_vecs(
		distmx, L, M, m,
		pens, pensids, mens, mensids);

//void chaq::fill_nen_ren_vecs(
//	cp_uint16_t pens,
//	cp_uint16_t mens,
//	cp_sid_t pensids,
//	cp_sid_t mensids,
//	uint L,
//	p_uint16_t nens,
//	p_uint16_t rens,
//	p_uint16_t nensids,
//	p_uint16_t rensids)
	chaq::fill_nen_ren_vecs(
		pens, mens, pensids, mensids, L,
		nens, rens, nensids, rensids);

	// Only "float" features, not aa, ss3 etc.
	const size_t bytes = L*sizeof(uint16_t);
	switch (fan)
		{
	case FAN_nendist:	memcpy(values, nensids, bytes); break;
	case FAN_rendist:	memcpy(values, rensids, bytes); break;
	case FAN_pendist:	memcpy(values, pensids, bytes); break;
	case FAN_mendist:	memcpy(values, mensids, bytes); break;

	case FAN_pm:
		{
		// Special-case hack, convert 8- to 16-bit.
		uint8_t *codeseq = myalloc(uint8_t, L);
		chaq::get_pm_codeseq(pensids, mensids, L, codeseq);
		for (uint i = 0; i < L; ++i)
			values[i] = codeseq[i];
		myfree(codeseq);
		break;
		}

	default:	Die("update_counts(%s)", FAN2str(fan));
		}

	myfree(distmx);
	myfree(pens);
	myfree(mens);
	myfree(nens);
	myfree(rens);
	myfree(pensids);
	myfree(mensids);
	myfree(nensids);
	myfree(rensids);
	}

static void update_counts(
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	uint16_t *counts)
	{
	const uint L = chain->get_length();
	uint16_t *values = myalloc(uint16_t, L);
	for (uint i = 0; i < L; ++i)
		values[i] = UINT16_MAX-1;
	get_values(chain, fan, alpha_size, values);

	for (uint i = 0; i < L; ++i)
		{
		uint16_t value = values[i];
		asserta(value != UINT16_MAX-1);
		counts[value] += 1;
		}

	myfree(values);
	}

static void make_charseq(
	const flat_chain_t *chain,
	FAN fan,
	uint8_t alpha_size,
	const uint16_t *thresholds,
	uint16_t undef_value,
	char *charseq)
	{
	const uint L = chain->get_length();

	uint16_t *values = myalloc(uint16_t, L);
	get_values(chain, fan, alpha_size, values);
	for (uint i = 0; i < L; ++i)
		{
		uint16_t value = values[i];
		if (value == UINT16_MAX)
			value = undef_value;
		uint8_t code = get_bin(value, alpha_size, thresholds);
		assert(code < alpha_size);
		charseq[i] = g_LetterToCharMu[code];
		}
	myfree(values);
	}

void cmd_flat_pm()
	{
	asserta(optset_fasta);
	const string &chainfn = g_Arg1;
	vector<vector<uint8_t> > codeseqs;
	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);
	FILE *f = CreateStdioFile(opt(fasta));
	for (uint i = 0; i < nchain; ++i)
		{
		uint L = chains[i]->get_length();
		uint16_t *values = myalloc(uint16_t, L);
		get_values(chains[i], FAN_pm, 2, values);
		string Seq;
		Seq.resize(L);
		for (uint i = 0; i < L; ++i)
			{
			uint v = values[i];
			asserta(v == 0 || v == 1);
			Seq[i] = (v == 0 ? 'A' : 'B');
			}
		SeqToFasta(f, chains[i]->m_label, Seq);
		myfree(values);
		}
	CloseStdioFile(f);
	}

void cmd_flat_aan()
	{
	asserta(optset_fasta);
	asserta(optset_alpha_size);
	const uint alpha_size = opt(alpha_size);
	const string &chainfn = g_Arg1;
	vector<vector<uint8_t> > codeseqs;
	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);
	FILE *f = CreateStdioFile(opt(fasta));
	for (uint i = 0; i < nchain; ++i)
		{
		const flat_chain_t &chain = *chains[i];
		uint L = chain.get_length();
		const char *aacharseq = chain.m_aa->m_data;
		uint8_t *codeseq = myalloc(uint8_t, L);
		switch (alpha_size)
			{
		case 3:
			chaq::get_aa3_codeseq(aacharseq, L, codeseq);
			break;

		case 4:
			chaq::get_aa4_codeseq(aacharseq, L, codeseq);
			break;

		default:
			Die("aan alpha_size %u", alpha_size);
			}

		string Seq;
		for (uint i = 0; i < L; ++i)
			{
			uint v = codeseq[i];
			Seq.push_back(g_LetterToCharMu[v]);
			}
		SeqToFasta(f, chains[i]->m_label, Seq);
		myfree(codeseq);
		}
	CloseStdioFile(f);
	}

void cmd_flat_quantize()
	{
	const string &chainfn = g_Arg1;
	vector<vector<uint8_t> > codeseqs;
	asserta(optset_alpha_size);
	const uint alpha_size = opt(alpha_size);

	asserta(optset_feature);
	const string feature_name = opt(feature);
	FAN fan = str2FAN(feature_name.c_str());

	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);

	vector<uint16_t> counts(UINT16_MAX+1);
	uint nbad = 0;
	for (uint i = 0; i < nchain; ++i)
		update_counts(chains[i], fan, alpha_size, counts.data());

	const QuantizeResult QR =
		quantize_histogram_equal_mass_dp(counts, alpha_size);

	double ideal_bin_size = double(QR.sum_count)/alpha_size;
	double sum_abs_diff = 0;
	const vector<uint16_t> &ts = QR.thresholds;

	string cmd;
	GetCmdLine(cmd);
	time_t t = time(0);
	char timeString[16];
	strftime(timeString, size(timeString), "%Y-%m-%d", gmtime(&t));

	ProgressLog("X  %7.7s  %7.7s  %7.7s\n", "Thresh.", "Size", "Diff");
	for (uint i = 0; i < alpha_size; ++i)
		{
		uint64_t bin_size = QR.bin_counts[i];
		double diff = double(bin_size) - ideal_bin_size;
		sum_abs_diff += abs(diff);
		uint16_t t = (i + 1 == alpha_size ? 0 : ts[i]);
		ProgressLog("%c", g_LetterToCharMu[i]);
		if (i + 1 < alpha_size)
			ProgressLog("  %7u", ts[i]);
		else
			ProgressLog("  %7.7s", "");
		ProgressLog("  %7u  %7.0f", bin_size, diff);
		ProgressLog("\n");
		}

	double mean_diff = sum_abs_diff/alpha_size;
	ProgressLog("Ideal %.1f, mean diff %.1f (%.1f%%) median %u\n",
		ideal_bin_size,
		mean_diff,
		GetPct(mean_diff, ideal_bin_size),
		QR.median_value);

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt(output));
		fprintf(f, "# %s\n", cmd.c_str());
		fprintf(f, "# [%s] %s\n", GIT_HASH, timeString);
		fprintf(f, "bins\t%u\n", alpha_size);
		for (uint i = 0; i + 1 < alpha_size; ++i)
			fprintf(f, "%u\t%u\n", i, ts[i]);
		CloseStdioFile(f);
		}

	if (optset_fasta)
		{
		FILE *f = CreateStdioFile(opt(fasta));
		for (uint i = 0; i < nchain; ++i)
			{
			uint L = chains[i]->get_length();
			string Seq;
			Seq.resize(L);
			make_charseq(chains[i], fan, alpha_size,
				ts.data(), QR.median_value, Seq.data());
			SeqToFasta(f, chains[i]->m_label, Seq);
			}
		CloseStdioFile(f);
		}
	}

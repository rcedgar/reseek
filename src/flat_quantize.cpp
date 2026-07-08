#include "myutils.h"
#include "chaq.h"
#include "flat_chain.h"
#include "flat_params.h"
#include "alpha.h"
#include "quantize.h"

void make_alpha_fnprefix(FAN fan, uint alpha_size, string &fnprefix)
	{
	string alphadir = optset_alphadir ? opt(alphadir) : "../alphadir";
	Dirize(alphadir);
	Ps(fnprefix, "%s%s%u", alphadir.c_str(), FAN2str(fan), alpha_size);
	}

static void update_counts(
	const flat_chain_t *chain,
	FAN fan,
	uint alpha_size,
	uint16_t *counts)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	const uint L = chain->get_length();
	uint16_t *values = myalloc(uint16_t, L);
	for (uint i = 0; i < L; ++i)
		values[i] = UINT16_MAX;
	chaq::slow_get_values(chain, fan, alpha_size, values);

	for (uint i = 0; i < L; ++i)
		{
		uint16_t value = values[i];
		counts[value] += 1;
		}

	myfree(values);
	}

void cmd_flat_quantize()
	{
	asserta(optset_alphadir);

	asserta(!optset_output);
	asserta(!optset_output2);
	asserta(!optset_fasta);

	flat_params params;

	uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	const string &chainfn = g_Arg1;
	vector<vector<uint8_t> > codeseqs;
	asserta(optset_alpha_size);
	const uint alpha_size = opt(alpha_size);
	if (optset_bandwidth) M = opt(bandwidth);

	asserta(optset_feature);
	const string feature_name = opt(feature);
	FAN fan = str2FAN(feature_name.c_str());

	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);

	vector<uint16_t> counts(UINT16_MAX+1);
	uint nbad = 0;
	for (uint i = 0; i < nchain; ++i)
		{
		ProgressStep(i, nchain, "Counting");
		update_counts(chains[i], fan, alpha_size, counts.data());
		}

	ProgressLog("Quantize...");
	const QuantizeResult QR =
		quantize_histogram_equal_mass_dp(counts, alpha_size);
	ProgressLog(" done.\n");

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

	// quantize
	string fnprefix;
	make_alpha_fnprefix(fan, alpha_size, fnprefix);

	const string fnq = fnprefix + ".quantize";
	ProgressLog("quantize: %s\n", fnq.c_str());
	FILE *fq = CreateStdioFile(fnq);
	fprintf(fq, "# %s\n", cmd.c_str());
	fprintf(fq, "# [%s] %s\n", GIT_HASH, timeString);
	fprintf(fq, "fan\t%s\n", FAN2str(fan));
	fprintf(fq, "alpha_size\t%u\n", alpha_size);
	fprintf(fq, "median\t%u\n", QR.median_value);
	fprintf(fq, "thresholds\t%u", alpha_size-1);
	for (uint i = 0; i + 1 < alpha_size; ++i)
		fprintf(fq, "\t%u", ts[i]);
	fprintf(fq, "\n");
	CloseStdioFile(fq);
	fq = 0;

	const string fna = fnprefix + ".fasta";
	ProgressLog("fasta %s\n", fna.c_str());
	FILE *ffa = CreateStdioFile(fna);
	for (uint i = 0; i < nchain; ++i)
		{
		uint L = chains[i]->get_length();
		char *Seq = myalloc(char, L);
		chaq::slow_get_charseq_binned(params, chains[i], fan, alpha_size,
			ts.data(), QR.median_value, Seq);
		SeqToFasta(ffa, chains[i]->m_label.c_str(), Seq, L);
		myfree(Seq);
		}
	CloseStdioFile(ffa);
	ffa = 0;
	}

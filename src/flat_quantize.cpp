#include "myutils.h"
#include "chaq.h"
#include "flat_chain.h"
#include "alpha.h"
#include "quantize.h"

static void update_counts(
	const flat_chain_t *chain,
	const string &feature,
	uint16_t *counts)
	{
	const uint L = chain->get_length();

	const uint M = 48;
	const uint m = 12;
	chaindistmx_t *dm;
	chaq::create_distmx(chain, dm, M);

	nnvec_t *nnvec;
	sidvec_t *nndisvec;
	chaq::create_nenvec(dm->m_data, M, L, m, nnvec, nndisvec);
	const sid_t *nns = nndisvec->m_data;
	for (uint i = 0; i < L; ++i)
		counts[nns[i]] += 1;
	}

static void make_charseq(
	const flat_chain_t *chain,
	const string &feature,
	uint8_t alpha_size,
	const uint16_t *thresholds,
	char *codeseq)
	{
	const uint L = chain->get_length();

	const uint M = 48;
	const uint m = 12;
	chaindistmx_t *dm;
	chaq::create_distmx(chain, dm, M);

	nnvec_t *nnvec;
	sidvec_t *nndistvec;
	chaq::create_nenvec(dm->m_data, M, L, m, nnvec, nndistvec);
	const sid_t *nns = nndistvec->m_data;
	for (uint i = 0; i < L; ++i)
		{
		uint16_t nndist = nns[i];
		uint8_t code = get_bin(nndist, alpha_size, thresholds);
		assert(code < alpha_size);
		codeseq[i] = g_LetterToCharMu[code];
		}
	}

void cmd_flat_quantize()
	{
	const string &chainfn = g_Arg1;
	vector<vector<uint8_t> > codeseqs;
	asserta(optset_alpha_size);
	const uint alpha_size = opt(alpha_size);

	asserta(optset_feature);
	const string feature = opt(feature);

	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);

	vector<uint16_t> counts(UINT16_MAX+1);
	uint nbad = 0;
	uint N = 0;
	for (uint i = 0; i < nchain; ++i)
		{
		update_counts(chains[i], feature, counts.data());
		N += chains[i]->get_length();
		}

	QuantizeResult QR = quantize_histogram_equal_mass_dp(counts, alpha_size);

	double ideal_bin_size = double(N)/alpha_size;
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
	ProgressLog("Ideal %.1f, mean diff %.1f  (%.1f%%)\n",
		ideal_bin_size, mean_diff, GetPct(mean_diff, ideal_bin_size));

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
			make_charseq(chains[i], feature, alpha_size, ts.data(), Seq.data());
			SeqToFasta(f, chains[i]->m_label, Seq);
			}
		CloseStdioFile(f);
		}
	}
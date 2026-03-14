#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"

void trunc_label(const string &Label,
	string &TruncatedLabel)
	{
	TruncatedLabel = Label;
	size_t n = TruncatedLabel.find(' ');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	n = TruncatedLabel.find('|');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	n = TruncatedLabel.find('/');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	}

void get_freqs_from_unaln(const string &fafn, uint alpha_size,
	vector<double> &freqs)
	{
	freqs.clear();
	freqs.resize(alpha_size);
	SeqDB db_fa;
	db_fa.FromFasta(fafn, false);
	const uint8_t *char2letter = (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);
	const uint nseq = db_fa.GetSeqCount();
	vector<uint> counts(alpha_size);
	uint N = 0;
	for (uint i = 0; i < nseq; ++i)
		{
		const uint L = db_fa.GetSeqLength(i);
		const string &seq = db_fa.GetSeq(i);
		for (uint j = 0; j < L; ++j)
			{
			uint8_t code = char2letter[seq[j]];
			if (code < alpha_size)
				{
				++counts[code];
				++N;
				}
			}
		}
	double sumfreq = 0;
	for (uint i = 0; i < alpha_size; ++i)
		{
		double freq = counts[i]/double(N);
		sumfreq += freq;
		freqs[i] = freq;
		}
	asserta(sumfreq > 0.99 && sumfreq < 1.01);

	Log("\nfreqs (unaln)\n");
	for (uint i = 0; i < alpha_size; ++i)
		Log("[%2u]  %6.4f\n", i, freqs[i]);
	}

void read_feature_fa_and_fa2(
	const string &fafn,
	const string &fa2fn,
	uint min_length,
	uint alpha_size,
	vector<uint8_t> &code1s,
	vector<uint8_t> &code2s)
	{
	code1s.clear();
	code2s.clear();

	const uint8_t *char2letter = (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);

	SeqDB db_fa, db_fa2;
	db_fa.FromFasta(fafn, false);
	db_fa2.FromFasta(fa2fn, true);

	const uint nfa = db_fa.GetSeqCount();
	const uint nfa2 = db_fa2.GetSeqCount();
	asserta(nfa2%2 == 0);
	map<string, uint> label2featseqidx;
	for (uint i = 0; i < nfa; ++i)
		{
		const uint L = db_fa.GetSeqLength(i);
		if (L < min_length)
			continue;
		const string &full_label = db_fa.GetLabel(i);
		string label;
		trunc_label(full_label, label);
		const string &featseq = db_fa.GetSeq(i);
		label2featseqidx[label] = i;
		}

	code1s.reserve(nfa2*400);
	code2s.reserve(nfa2*400);
	uint length_mismatch_count = 0;
	uint npairs = 0;
	uint nmissing = 0;
	uint bad_letters = 0;
	for (uint i = 0; i < nfa2; i += 2)
		{
		uint ncols = db_fa2.GetSeqLength(i);
		asserta(db_fa2.GetSeqLength(i + 1) == ncols);
		const string &full_label1 = db_fa2.GetLabel(i);
		const string &full_label2 = db_fa2.GetLabel(i+1);
		string label1, label2;
		trunc_label(full_label1, label1);
		trunc_label(full_label2, label2);
		map<string, uint>::const_iterator iter1 = label2featseqidx.find(label1);
		map<string, uint>::const_iterator iter2 = label2featseqidx.find(label2);
		if (iter1 == label2featseqidx.end() || iter2 == label2featseqidx.end())
			{
			++nmissing;
			continue;
			}
		uint featseqidx1 = iter1->second;
		uint featseqidx2 = iter2->second;
		uint L1 = db_fa.GetSeqLength(featseqidx1);
		uint L2 = db_fa.GetSeqLength(featseqidx2);
		const string &row1 = db_fa2.GetSeq(i);
		const string &row2 = db_fa2.GetSeq(i+1);
		asserta(SIZE(row1) == ncols);
		asserta(SIZE(row2) == ncols);
		uint pos1 = 0;
		uint pos2 = 0;
		for (uint colidx = 0; colidx < ncols; ++colidx)
			{
			char c1 = row1[colidx];
			char c2 = row2[colidx];
			if (!isgap(c1)) ++pos1;
			if (!isgap(c2)) ++pos2;
			}
		if (pos1 != L1 || pos2 != L2)
			{
			string useq1, useq2;
			db_fa2.GetUngappedSeq(i, useq1);
			db_fa2.GetUngappedSeq(i+1, useq2);
			Log("\nmismatch pos1=%u L1=%u, pos2=%u L2=%u\n", pos1, L1, pos2, L2);
			Log("%s %s\n", db_fa2.GetSeq(i).c_str(), label1.c_str());
			Log("%s %s\n", useq1.c_str(), label1.c_str());
			Log("%s %s\n", db_fa.GetSeq(featseqidx1).c_str(), label1.c_str());
			Log("\n");
			Log("%s %s\n", db_fa2.GetSeq(i+1).c_str(), label2.c_str());
			Log("%s %s\n", useq2.c_str(), label2.c_str());
			Log("%s %s\n", db_fa.GetSeq(featseqidx2).c_str(), label2.c_str());
			++length_mismatch_count;
			continue;
			}
		++npairs;
		pos1 = 0;
		pos2 = 0;
		const string &featseq1 = db_fa.GetSeq(featseqidx1);
		const string &featseq2 = db_fa.GetSeq(featseqidx2);
		asserta(SIZE(featseq1) == L1);
		asserta(SIZE(featseq2) == L2);
		//string featrow1;//@@
		//string featrow2;//@@
		for (uint colidx = 0; colidx < ncols; ++colidx)
			{
			char rowc1 = row1[colidx];
			char rowc2 = row2[colidx];
			if (isupper(rowc1) && isupper(rowc2))
				{
				uint8_t featc1 = featseq1[pos1];
				uint8_t featc2 = featseq2[pos2];
				uint8_t code1 = char2letter[featc1];
				uint8_t code2 = char2letter[featc2];
				if (code1 < alpha_size && code2 < alpha_size)
					{
					//featrow1 += featc1;
					//featrow2 += featc2;
					code1s.push_back(code1);
					code2s.push_back(code2);
					}
				else
					++bad_letters;
				}
			if (!isgap(rowc1)) ++pos1;
			if (!isgap(rowc2)) ++pos2;
			}
		asserta(pos1 == L1);
		asserta(pos2 == L2);
		//Log("\n");
		//Log("%s, %s\n", label1.c_str(), label2.c_str());
		//Log("%s\n", row1.c_str());
		//Log("%s\n", row2.c_str());
		//Log("%s\n", featrow1.c_str());
		//Log("%s\n", featrow2.c_str());
		}
	ProgressLog("%u seq pairs, %s letter pairs, %u length mismatches, %u missing, %u bad\n",
		npairs, IntToStr(SIZE(code1s)), length_mismatch_count, nmissing, bad_letters);
	}

void get_countmx_from_code_pairs(
	const vector<uint8_t> &code1s,
	const vector<uint8_t> &code2s,
	uint alpha_size,
	vector<vector<uint> > &countmx)
	{
	const uint n = SIZE(code1s);
	asserta(SIZE(code2s) == n);
	countmx.clear();
	countmx.resize(alpha_size);
	for (uint i = 0; i < alpha_size; ++i)
		countmx[i].resize(alpha_size);
	uint M = 0;
	for (uint i = 0; i < n; ++i)
		{
		uint8_t code1 = code1s[i];
		uint8_t code2 = code2s[i];
		asserta(code1 < alpha_size && code2 < alpha_size);
		countmx[code1][code2] += 1;
		countmx[code2][code1] += 1;
		M += 1;
		}
	asserta(M == n);

	uint N = 0;
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = 0; j < alpha_size; ++j)
			{
			uint n = countmx[i][j];
			Log("  %7u", n);
			N += n;
			}
		Log("\n");
		}
	asserta(N == 2*n);
	}

void get_marginal_freqs_from_code_countsmx(
	const vector<vector<uint> > &countmx,
	vector<double> &freqs)
	{
	freqs.clear();
	const uint alpha_size = SIZE(countmx);
	vector<uint> ns;
	uint N = 0;
	for (uint i = 0; i < alpha_size; ++i)
		{
		uint n = 0;
		for (uint j = 0; j < alpha_size; ++j)
			{
			asserta(countmx[i][j] == countmx[j][i]);
			n += countmx[i][j];
			}
		N += n;
		ns.push_back(n);
		}
	double sumfreq = 0;
	for (uint i = 0; i < alpha_size; ++i)
		{
		double freq = double(ns[i])/N;
		freqs.push_back(freq);
		sumfreq += freq;
		}
	asserta(sumfreq > 0.99 && sumfreq < 1.01);

	Log("\nfreqs (cols)\n");
	for (uint i = 0; i < alpha_size; ++i)
		Log("[%2u]  %6.4f\n", i, freqs[i]);
	}

void get_joint_freqmx_from_code_countsmx(
	const vector<vector<uint> > &countmx,
	vector<vector<double> > &freqmx)
	{
	freqmx.clear();
	const uint alpha_size = SIZE(countmx);
	freqmx.resize(alpha_size);
	for (uint i = 0; i < alpha_size; ++i)
		freqmx[i].resize(alpha_size);

	uint N = 0;
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = 0; j < alpha_size; ++j)
			{
			asserta(countmx[i][j] == countmx[j][i]);
			N += countmx[i][j];
			}
		}

	double sumfreq = 0;
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = 0; j < alpha_size; ++j)
			{
			asserta(countmx[i][j] == countmx[j][i]);
			uint n = countmx[i][j];
			double freq = double(n)/N;
			freqmx[i][j] = freq;
			sumfreq += freq;
			}
		}
	asserta(sumfreq > 0.99 && sumfreq < 1.01);

	Log("\nfreqmx\n");
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = 0; j < alpha_size; ++j)
			Log("  %8.2g", freqmx[i][j]);
		Log("\n");
		}
	}

double log_base(double x, double base)
	{
	return log(x)/log(base);
	}

void get_logoddsmx_from_freqs(
	const vector<double> &freqs,
	const vector<vector<double> > &freqmx,
	vector<vector<double> > &logoddsmx,
	const string &units)
	{
	uint alpha_size = SIZE(freqs);
	asserta(SIZE(freqs) == alpha_size);
	asserta(SIZE(freqmx) == alpha_size);
	logoddsmx.resize(alpha_size);
	for (uint i = 0; i < alpha_size; ++i)
		logoddsmx[i].resize(alpha_size);

	for (uint i = 0; i < alpha_size; ++i)
		{
		double f_i = freqs[i];
		for (uint j = 0; j < alpha_size; ++j)
			{
			double f_j = freqs[j];
			double f_ij = freqmx[i][j];
			asserta(f_ij > 1e-6);
			asserta(feq(freqmx[j][i], f_ij));
			if (units == "bits")
				logoddsmx[i][j] = log2(f_ij) - log2(f_i) - log2(f_j);
			else if (units == "nats")
				logoddsmx[i][j] = log(f_ij) - log(f_i) - log(f_j);
			else
				Die("units=%s", units.c_str());
			}
		}
	Log("\nlogodds\n");
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = 0; j < alpha_size; ++j)
			Log("  %8.2g", logoddsmx[i][j]);
		Log("\n");
		}
	}

void scale_logoddsmx(vector<vector<double> > &logoddsmx, uint scale)
	{
	const uint alpha_size = SIZE(logoddsmx);
	for (uint i = 0; i < alpha_size; ++i)
		for (uint j = 0; j < alpha_size; ++j)
			logoddsmx[i][j] = round(scale*logoddsmx[i][j]);
	}

void write_logoddsmx(FILE *f,
	const vector<vector<double> > &logoddsmx,
	bool asintegers)
	{
	if (f == 0)
		return;
	string cmd;
	GetCmdLine(cmd);
	const uint alpha_size = SIZE(logoddsmx);
	time_t t = time(0);
	char timeString[16];
	strftime(timeString, size(timeString), "%Y-%m-%d", gmtime(&t));
	fprintf(f, "# %s\n", cmd.c_str());
	fprintf(f, "# [%s] %s\n", GIT_HASH, timeString);
	fprintf(f, "logodds\t%u\n", alpha_size);
	if (asintegers)
		{
		for (uint i = 0; i < alpha_size; ++i)
			{
			fprintf(f, "%u", i);
			for (uint j = 0; j < alpha_size; ++j)
				fprintf(f, "\t%d", int(round(logoddsmx[i][j])));
			fprintf(f, "\n");
			}
		}
	else
		{
		for (uint i = 0; i < alpha_size; ++i)
			{
			fprintf(f, "%u", i);
			for (uint j = 0; j < alpha_size; ++j)
				fprintf(f, "\t%.4g", logoddsmx[i][j]);
			}
		fprintf(f, "\n");
		}
	}

void cmd_flat_train_discrete()
	{
	const string &fafn = g_Arg1;
	const string &fa2fn = opt(fasta2_tp);
	vector<vector<uint8_t> > codeseqs;
	const uint min_length = 100;
	asserta(optset_alpha_size);
	const uint alpha_size = opt(alpha_size);

	vector<uint8_t> code1s;
	vector<uint8_t> code2s;
	read_feature_fa_and_fa2(fafn, fa2fn, min_length, alpha_size,
		code1s, code2s);

	vector<vector<uint> > countmx;
	get_countmx_from_code_pairs(code1s, code2s, alpha_size, countmx);

	vector<vector<double> > freqmx;
	asserta(optset_background_style);
	const string bs = opt(background_style);

	get_joint_freqmx_from_code_countsmx(countmx, freqmx);

	vector<double> freqs_cols;
	vector<double> freqs_unaln;
	get_marginal_freqs_from_code_countsmx(countmx, freqs_cols);
	get_freqs_from_unaln(fafn, alpha_size, freqs_unaln);

	for (uint i = 0; i < alpha_size; ++i)
		{
		double fc = freqs_cols[i];
		double fu = freqs_unaln[i];
		double r = fc/fu;
		Log("%6.4f  %6.4f  %6.4f\n", fc, fu, r);
		}

	vector<vector<double> > logoddsmx;
	string units("bits");
	if (optset_units)
		units = opt(units);

	if (bs == "cols")
		get_logoddsmx_from_freqs(freqs_cols, freqmx, logoddsmx, units);
	else if (bs == "unaln")
		get_logoddsmx_from_freqs(freqs_unaln, freqmx, logoddsmx, units);
	else
		Die("-background_style %s", bs.c_str());

	uint scale = 1;
	if (optset_scale)
		scale_logoddsmx(logoddsmx, opt(scale));
	FILE *f = CreateStdioFile(opt(output));
	bool asintegers = optset_scale;
	write_logoddsmx(f, logoddsmx, asintegers);
	CloseStdioFile(f);
	}

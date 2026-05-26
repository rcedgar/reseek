#include "myutils.h"
#include "flat_helpers.h"
#include "hexintseq.h"
#include "sort.h"
#include "seqdb.h"
#include "alpha.h"
#include <numeric>

static uint get_maxreducedcode(
	const vector<uint> &fullcode2reducedcode)
	{
	const uint fullAS = uint(fullcode2reducedcode.size());
	uint maxreducedcode = 0;
	for (uint full_code = 0; full_code < fullAS; ++full_code)
		{
		uint reduced_code = fullcode2reducedcode[full_code];
		maxreducedcode = max(reduced_code, maxreducedcode);
		}
	asserta(maxreducedcode > 0);
	return maxreducedcode;
	}

static uint get_reducedAS(
	const vector<uint> &fullcode2reducedcode)
	{
	return 1 + get_maxreducedcode(fullcode2reducedcode);
	}

static void invert_map(
	const vector<uint> &fullcode2reducedcode,
	vector<vector<uint> > &reducedcode2full_codes)
	{
	const uint fullAS = uint(fullcode2reducedcode.size());
	const uint reducedAS = get_reducedAS(fullcode2reducedcode);
	reducedcode2full_codes.clear();
	reducedcode2full_codes.resize(reducedAS);
	for (uint fullcode = 0; fullcode < fullAS; ++fullcode)
		{
		uint reducedcode = fullcode2reducedcode[fullcode];
		assert(reducedcode < reducedAS);
		reducedcode2full_codes[reducedcode].push_back(fullcode);
		}
	}

static void log_map(
	const vector<uint> &fullcode2reducedcode,
	double H)
	{
	uint fullAS = uint(fullcode2reducedcode.size());
	uint reducedAS = get_reducedAS(fullcode2reducedcode);
	vector<vector<uint> > reducedcode2full_codes;
	invert_map(fullcode2reducedcode, reducedcode2full_codes);
	ProgressLog("[%3u] ", reducedAS);
	ProgressLog(" %.4f", H);
	Log(" | ");
	for (uint i = 0; i < reducedAS; ++i)
		{
		const vector<uint> &fullcodes =
			reducedcode2full_codes[i];
		Log(" (");
		for (auto fullcode : fullcodes)
			{
			Log(" %u", fullcode);
			}
		Log(" )");
		}
	ProgressLog("\n");
	}

// After merge, merged pair -> newAS-1; 
//   other old reduced codes -> 0..newAS-2.
static void merge_old2new(
	uint reducedAS, uint i, uint j, vector<uint> &old2new)
	{
	asserta(i < reducedAS && j < reducedAS);
	asserta(i != j);
	const uint newAS = reducedAS - 1;
	old2new.assign(reducedAS, UINT_MAX);
	uint idx = 0;
	for (uint old = 0; old < reducedAS; ++old)
		{
		if (old == i || old == j)
			continue;
		old2new[old] = idx++;
		}
	asserta(idx + 1 == newAS);
	old2new[i] = newAS - 1;
	old2new[j] = newAS - 1;
	}

static double getH(
	const vector<double> &fullfreqmx,
	const vector<uint> &fullcode2reducedcode)
	{
	const uint reducedAS = get_reducedAS(fullcode2reducedcode);
	const uint fullAS = uint(fullcode2reducedcode.size());
	asserta(fullfreqmx.size() == fullAS*fullAS);

	vector<double> reducedfreqmx(reducedAS*reducedAS);
	for (uint fullcode1 = 0; fullcode1 < fullAS; ++fullcode1)
		{
		uint r1 = fullcode2reducedcode[fullcode1];
		for (uint fullcode2 = 0; fullcode2 < fullAS; ++fullcode2)
			{
			uint r2 = fullcode2reducedcode[fullcode2];
			double freq = fullfreqmx[fullcode1*fullAS + fullcode2];
			reducedfreqmx[r1*reducedAS + r2] += freq;
			}
		}

	vector<double> logoddsmx;
	get_logoddsmx_from_flat_freqmx(
		reducedfreqmx, reducedAS, logoddsmx);
	double H = get_relative_entropy_flat(
		reducedfreqmx, logoddsmx, reducedAS);

	return H;
	}

static double merge_letters(
	const vector<double> &fullfreqmx,
	const vector<uint> &fullcode2reducedcode,
	uint i, uint j)
	{
	const uint maxreducedcode = get_maxreducedcode(fullcode2reducedcode);
	const uint reducedAS = maxreducedcode + 1;
	const uint fullAS = uint(fullcode2reducedcode.size());
	asserta(i < reducedAS && j < reducedAS);
	asserta(i != j);
	asserta(fullfreqmx.size() == fullAS*fullAS);

	const uint newAS = reducedAS - 1;
	vector<uint> old2new;
	merge_old2new(reducedAS, i, j, old2new);

	vector<double> reducedfreqmx(newAS*newAS);

	for (uint fullcode1 = 0; fullcode1 < fullAS; ++fullcode1)
		{
		uint r1 = old2new[fullcode2reducedcode[fullcode1]];
		for (uint fullcode2 = 0; fullcode2 < fullAS; ++fullcode2)
			{
			uint r2 = old2new[fullcode2reducedcode[fullcode2]];
			double freq = fullfreqmx[fullcode1*fullAS + fullcode2];
			reducedfreqmx[r1*newAS + r2] += freq;
			}
		}

	vector<double> logoddsmx;
	get_logoddsmx_from_flat_freqmx(
		reducedfreqmx, newAS, logoddsmx);
	double H = get_relative_entropy_flat(
		reducedfreqmx, logoddsmx, newAS);

	return H;
	}

static double find_best_ij(
	const vector<double> &in_freqmx,
	uint alpha_size,
	uint &best_i,
	uint &best_j,
	const vector<uint> &fullcode2reducedcode)
	{
	asserta(alpha_size >= 2);
	bool have = false;
	double best_H = 0;
	best_i = UINT_MAX;
	best_j = UINT_MAX;
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = i+1; j < alpha_size; ++j)
			{
			double H = merge_letters(
				in_freqmx, fullcode2reducedcode, i, j);
			if (!have || H > best_H)
				{
				have = true;
				best_H = H;
				best_i = i;
				best_j = j;
				}
			}
		}
	return best_H;
	}

static void validate_map(const vector<uint> &fullcode2reducedcode)
	{
	uint fullAS = uint(fullcode2reducedcode.size());
	uint maxreducedcode = get_maxreducedcode(fullcode2reducedcode);
	uint reducedAS = maxreducedcode + 1;
	vector<bool> found(reducedAS);
	uint nfound = 0;
	for (uint full_code = 0; full_code < fullAS ; ++full_code)
		{
		uint reduced_code = fullcode2reducedcode[full_code];
		if (!found[reduced_code])
			{
			found[reduced_code] = true;
			++nfound;
			}
		}
	asserta(nfound == reducedAS);
	}

static void upd_map(
	vector<uint> &fullcode2reducedcode,
	uint best_i, uint best_j)
	{
	validate_map(fullcode2reducedcode);

	uint fullAS = uint(fullcode2reducedcode.size());
	uint maxreducedcode = get_maxreducedcode(fullcode2reducedcode);
	const uint reducedAS = maxreducedcode + 1;

	asserta(best_i < reducedAS && best_j < reducedAS);
	asserta(best_i != best_j);

	const uint newAS = reducedAS - 1;
	vector<uint> old2new;
	merge_old2new(reducedAS, best_i, best_j, old2new);

	for (uint fullcode = 0; fullcode < fullAS; ++fullcode)
		{
		uint oldcode = fullcode2reducedcode[fullcode];
		asserta(oldcode < reducedAS);
		fullcode2reducedcode[fullcode] = old2new[oldcode];
		asserta(fullcode2reducedcode[fullcode] < newAS);
		}

	validate_map(fullcode2reducedcode);
	}

void peturb(
	const vector<uint> &fullcode2reducedcode,
	vector<uint> &fc)
	{
	const uint fullAS = uint(fullcode2reducedcode.size());
	const uint reducedAS = get_reducedAS(fullcode2reducedcode);
	uint r = randu32()%100;
	uint nmut = 1;
	if (r > 80)
		nmut = 2;
	if (r > 95)
		nmut = 3;

	fc = fullcode2reducedcode;
	for (uint mutidx = 0; mutidx < nmut; ++mutidx)
		{
		switch (randu32()%2)
			{
		case 0: // swap
			{
			uint i = randu32()%reducedAS;
			uint j = randu32()%(reducedAS-1);
			if (i == j)
				++j;
			swap(fc[i], fc[j]);
			break;
			}

		case 1: // poke
			{
			uint reduced_code = randu32()%reducedAS;
			uint full_code = randu32()%fullAS;
			fc[full_code] = reduced_code;
			break;
			}
		
		default: asserta(false);
			}
		}
	}

static vector<vector<uint> > s_donevec;

static void add_done(const vector<uint> &fc)
	{
	s_donevec.push_back(fc);
	}

static bool check_isdone(const vector<uint> &fc)
	{
	for (auto v : s_donevec)
		if (fc == v)
			return true;
	return false;
	}

static bool normalize_fc(vector<uint> &fc)
	{
	uint reducedAS = get_reducedAS(fc);
	vector<uint> old2new(reducedAS, UINT_MAX);
	uint newcode = 0;
	for (uint i = 0; i < fc.size(); ++i)
		{
		uint reduced_code = fc[i];
		if (old2new[reduced_code] == UINT_MAX)
			old2new[reduced_code] = newcode++;
		}
	asserta(newcode <= reducedAS);
	if (newcode < reducedAS)
		return false;
	for (uint i = 0; i < fc.size(); ++i)
		{
		uint reduced_code = fc[i];
		fc[i] = old2new[reduced_code];
		}
	return true;
	}

void cmd_reduce_alphabet()
	{
	asserta(optset_alpha_size);
	asserta(!optset_n);

	const string &logoddsfn = g_Arg1;
	const uint reducedAS = opt(alpha_size);

	vector<vector<uint> > fcvec;
	vector<double> Hvec;

	vector<double> logoddsmx, in_freqmx;
	uint fullAS = UINT_MAX;
	if (logoddsfn[0] == '@')
		{
		string alpha_name = logoddsfn.substr(1);
		extern const vector<string> g_alpha_collect_lines;
		collect C;
		C.from_lines(g_alpha_collect_lines);
		string logoddsfn;
		Ps(logoddsfn, "%s.logodds", alpha_name.c_str());

		vector<float> logodds;
		const vector<string> &logodds_lines = C.get_lines(logoddsfn);
		fullAS = logodds_and_freqmx_from_lines(
			logodds_lines, logoddsmx, in_freqmx);
		}
	else
		fullAS = read_logodds_and_freqmx(logoddsfn, logoddsmx, in_freqmx);
	asserta(fullAS <= 36);
	const uint AS2 = fullAS*fullAS;
	asserta(logoddsmx.size() == AS2);
	asserta(in_freqmx.size() == AS2);
	asserta(reducedAS < fullAS);
	asserta(fullAS <= 36);

	SeqDB DB;
	if (optset_input)
		DB.FromFasta(opt(input));

	vector<double> freqs(fullAS);
	double sumfreqs = 0;
	for (uint i = 0; i < fullAS; ++i)
		{
		double freq = 0;
		for (uint j = 0; j < fullAS; ++j)
			freq += in_freqmx[fullAS*i + j];
		freqs[i] = freq;
		sumfreqs += freq;
		}
	asserta(sumfreqs > 0.99 && sumfreqs < 1.01);

	vector<uint> fullcode2reducedcode;
	for (uint i = 0; i < fullAS; ++i)
		fullcode2reducedcode.push_back(i);

	validate_map(fullcode2reducedcode);
	double H = 0;
	for (uint AS = fullAS; AS > reducedAS; --AS)
		{
		uint best_i, best_j;
		H = find_best_ij(in_freqmx,
			AS, best_i, best_j, fullcode2reducedcode);
		ProgressLog("[%3u]  %3u  %3u  H=%.4f\n",
			AS-1, best_i, best_j, H);

		upd_map(fullcode2reducedcode, best_i, best_j);
		}
	log_map(fullcode2reducedcode, H);

	double bestH = H;
	bool ok = normalize_fc(fullcode2reducedcode);
	asserta(ok);
	double H2 = getH(in_freqmx, fullcode2reducedcode);
	asserta(feq(H2, bestH));

	fcvec.push_back(fullcode2reducedcode);
	Hvec.push_back(bestH);

	const uint ITERS = 10000;
	vector<uint> fc;
	for (uint iter = 0; iter < ITERS; ++iter)
		{
		ProgressStep(iter, ITERS, "Perturbing");
		peturb(fullcode2reducedcode, fc);
		double H2 = getH(in_freqmx, fc);
		normalize_fc(fc);
		double H3 = getH(in_freqmx, fc);
		asserta(feq(H2, H3));

		fcvec.push_back(fc);
		Hvec.push_back(H2);

		if (H2 > bestH)
			{
			ProgressLog("|%5u| %6.4f <<<\n", iter, H2);
			fullcode2reducedcode = fc;
			bestH = H2;
			}
		}
	log_map(fullcode2reducedcode, H);

	for (uint i = 0; i < fullAS; ++i)
		fc[i] = randu32()%reducedAS;

	bestH = getH(in_freqmx, fc);
	log_map(fc, bestH);

	Progress("\n");
	Progress("\n");
	Progress("\n");
	const uint ITERS_SHUFFLE = 10000;
	for (uint iter = 0; iter < ITERS_SHUFFLE; ++iter)
		{
		ProgressStep(iter, ITERS_SHUFFLE, "Shuffle climb");
		peturb(fullcode2reducedcode, fc);
		double H2 = getH(in_freqmx, fc);

		fcvec.push_back(fc);
		Hvec.push_back(H2);

		if (H2 > bestH)
			{
			ProgressLog("|%5u| %6.4f <<<\n", iter, H2);
			fullcode2reducedcode = fc;
			bestH = H2;
			}
		}
	log_map(fullcode2reducedcode, H);

	const size_t n = fcvec.size();
	asserta(Hvec.size() == n);
	vector<uint> order(n);
#if 1
	iota(order.begin(), order.end(), 0u);
	sort(order.begin(), order.end(),
		[&](unsigned a, unsigned b) {
			return Hvec[a] > Hvec[b];
		});
#else
	QuickSortOrderDesc(Hvec.data(), uint(n), order.data());
#endif
	uint topn = 10;
	if (optset_topn)
		topn = opt(topn);
	if (topn > n)
		topn = uint(n);

	vector<vector<uint> > output_fcs;
	vector<double> output_Hs;
	double lastH = DBL_MAX;
	for (size_t k = 0; k < n; ++k)
		{
		uint i = order[k];
		if (s_donevec.size() >= topn)
			break;
		double H = Hvec[i];
		asserta(H <= lastH);
		lastH = H;
		const vector<uint> &fc = fcvec[i];
		bool isdone = check_isdone(fc);
		if (!isdone)
			{
			output_fcs.push_back(fc);
			output_Hs.push_back(H);
			add_done(fc);
			Log("top[%3u] = %.4f\n", uint(s_donevec.size()), H);
			}
		}

	if (optset_output)
		{
		FILE *fout = CreateStdioFile(opt(output));

		asserta(fullAS <= 36);
		const byte *letter2char = (fullAS == 20 ? g_LetterToCharAmino : g_LetterToCharMu);
		for (size_t i = 0; i < output_fcs.size(); ++i)
			{
			double H = output_Hs[i];
			const vector<uint> &fc = output_fcs[i];
			fprintf(fout, "%.4f", H);
			add_done(fc);
			for (uint full_code = 0; full_code < fullAS; ++full_code)
				fprintf(fout, "\t%u", fc[full_code]);
			vector<vector<uint> > inv;
			invert_map(fc, inv);
			asserta(inv.size() == reducedAS);
			vector<string> reduced_strings;
			for (uint reduced_code = 0; reduced_code < reducedAS; ++reduced_code)
				{
				string s;
				const vector<uint> &v = inv[reduced_code];
				fprintf(fout, "\t(");
				for (uint k = 0; k < v.size(); ++k)
					{
					if (k > 0)
						fprintf(fout, ",");

					uint code = v[k];
					asserta(code < fullAS);
					s += letter2char[code];
					fprintf(fout, "%u", code);
					}
				reduced_strings.push_back(s);
				fprintf(fout, ")");
				}
			fprintf(fout, "\t");
			for (uint reduced_code = 0; reduced_code < reducedAS; ++reduced_code)
				{
				if (reduced_code > 0)
					fprintf(fout, "-");
				fprintf(fout, "%s", reduced_strings[reduced_code].c_str());
				}

			fprintf(fout, "\n");
			}
		ProgressLog("%u written\n", uint(output_fcs.size()));
		CloseStdioFile(fout);
		}

	if (optset_output2)
		{
		asserta(optset_input);
		const uint nseq = DB.GetSeqCount();
		DB.ToLetters(g_CharToLetterMu);
		const string prefix = opt(output2);
		for (size_t i = 0; i < output_fcs.size(); ++i)
			{
			string fafn;
			char c = 'A' + uint8_t(i);
			Ps(fafn, "%s%c%u.fa", prefix.c_str(), c, reducedAS);
			ProgressLog("%s\n", fafn.c_str());
			FILE *f = CreateStdioFile(fafn);
			double H = output_Hs[i];
			const vector<uint> &fc = output_fcs[i];
			for (uint seqidx = 0; seqidx < nseq; ++seqidx)
				{
				const byte *byteseq = DB.GetByteSeq(seqidx);
				uint L = DB.GetSeqLength(seqidx);
				const string &label = DB.GetLabel(seqidx);
				string outseq;
				for (uint pos = 0; pos < L; ++pos)
					{
					uint8_t full_code = byteseq[pos];
					asserta(full_code < fullAS);
					uint8_t reduced_code = fc[full_code];
					char c = g_LetterToCharMu[reduced_code];
					outseq += c;
					}
				SeqToFasta(f, label, outseq);
				}
			CloseStdioFile(f);
			}
		}
	}

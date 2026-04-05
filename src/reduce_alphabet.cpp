#include "myutils.h"
#include "flat_helpers.h"

static double merge_letters(
	uint8_t i, uint8_t j, uint in_alpha_size,
	const vector<double> &in_freqmx,
	vector<double> &out_freqmx)
	{
	asserta(i < in_alpha_size);
	asserta(j < in_alpha_size);
	asserta(i != j);
	asserta(in_freqmx.size() == in_alpha_size*in_alpha_size);

	const uint out_alpha_size = in_alpha_size - 1;
	out_freqmx.clear();
	out_freqmx.resize(out_alpha_size*out_alpha_size);

	vector<uint> in_code_to_out_code;
	uint idx = 0;
	for (uint k = 0; k < in_alpha_size; ++k)
		{
		if (k == i || k == j)
			in_code_to_out_code.push_back(out_alpha_size-1);
		else
			in_code_to_out_code.push_back(idx++);
		}
	asserta(idx+1 == out_alpha_size);

	for (uint in_code1 = 0; in_code1 < in_alpha_size; ++in_code1)
		{
		uint out_code1 = in_code_to_out_code[in_code1];
		for (uint in_code2 = 0; in_code2 < in_alpha_size; ++in_code2)
			{
			uint out_code2 = in_code_to_out_code[in_code2];
			double freq = in_freqmx[in_code1*in_alpha_size + in_code2];
			out_freqmx[out_code1*out_alpha_size + out_code2] += freq;
			}
		}

	vector<double> logoddsmx;
	get_logoddsmx_from_flat_freqmx(
		out_freqmx, out_alpha_size, logoddsmx);
	double H = get_relative_entropy_flat(
		out_freqmx, logoddsmx, out_alpha_size);

	return H;
	}

static double find_best_ij(
	const vector<double> &in_freqmx,
	uint alpha_size,
	uint &best_i,
	uint &best_j,
	vector<double> &out_freqmx)
	{
	double best_H = 0;
	best_i = UINT_MAX;
	best_j = UINT_MAX;
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = i+1; j < alpha_size; ++j)
			{
			double H = merge_letters(
				i, j, alpha_size, in_freqmx, out_freqmx);
			if (H > best_H)
				{
				best_H = H;
				best_i = i;
				best_j = j;
				}
			}
		}
	double H = merge_letters(
		best_i, best_j, alpha_size, in_freqmx, out_freqmx);
	asserta(H == best_H);
	return H;
	}

static void validate_map(const vector<uint> &fullcode2reducedcode)
	{
	uint fullAS = uint(fullcode2reducedcode.size());
	uint maxreducedcode = 0;
	for (uint full_code = 0; full_code < fullAS ; ++full_code)
		{
		uint reduced_code = fullcode2reducedcode[full_code];
		maxreducedcode = max(reduced_code, maxreducedcode);
		}
	asserta(maxreducedcode > 0);
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
	uint maxreducedcode = 0;
	for (uint fullcode = 0; fullcode < fullAS; ++fullcode)
		maxreducedcode = max(maxreducedcode, fullcode2reducedcode[fullcode]);
	const uint reducedAS = maxreducedcode + 1;

	asserta(best_i < reducedAS && best_j < reducedAS);
	asserta(best_i != best_j);

	// Match merge_letters: merged pair -> last code (newAS-1); others compact 0..newAS-2
	const uint newAS = reducedAS - 1;
	vector<uint> old2new(reducedAS);
	for (uint old = 0; old < reducedAS; ++old)
		old2new[old] = UINT_MAX;
	uint idx = 0;
	for (uint old = 0; old < reducedAS; ++old)
		{
		if (old == best_i || old == best_j)
			continue;
		old2new[old] = idx++;
		}
	asserta(idx + 1 == newAS);
	old2new[best_i] = newAS - 1;
	old2new[best_j] = newAS - 1;

	for (uint fullcode = 0; fullcode < fullAS; ++fullcode)
		{
		uint oldcode = fullcode2reducedcode[fullcode];
		asserta(oldcode < reducedAS);
		fullcode2reducedcode[fullcode] = old2new[oldcode];
		asserta(fullcode2reducedcode[fullcode] < newAS);
		}

	validate_map(fullcode2reducedcode);
	}

void cmd_reduce_alphabet()
	{
	const string &logoddsfn = g_Arg1;
	asserta(optset_alpha_size);
	asserta(!optset_n);
	const uint out_alpha_size = opt(alpha_size);

	vector<double> logoddsmx, in_freqmx;
	const uint in_alpha_size = 
		read_logodds_and_freqmx(logoddsfn, logoddsmx, in_freqmx);
	const uint AS2 = in_alpha_size*in_alpha_size;
	asserta(logoddsmx.size() == AS2);
	asserta(in_freqmx.size() == AS2);
	asserta(out_alpha_size < in_alpha_size);

	vector<double> freqs(in_alpha_size);
	double sumfreqs = 0;
	for (uint i = 0; i < in_alpha_size; ++i)
		{
		double freq = 0;
		for (uint j = 0; j < in_alpha_size; ++j)
			freq += in_freqmx[in_alpha_size*i + j];
		freqs[i] = freq;
		sumfreqs += freq;
		}
	asserta(sumfreqs > 0.99 && sumfreqs < 1.01);
	double H = get_relative_entropy_flat(
		in_freqmx, logoddsmx, in_alpha_size);

	ProgressLog("Input ES=%.3g\n", H);

	vector<double> out_freqmx;
	vector<uint> fullcode2reducedcode;

	for (uint i = 0; i < in_alpha_size; ++i)
		fullcode2reducedcode.push_back(i);

	validate_map(fullcode2reducedcode);
	for (uint AS = in_alpha_size; AS >= out_alpha_size; --AS)
		{
		uint best_i, best_j;
		asserta(in_freqmx.size() == AS*AS);
		double H = find_best_ij(in_freqmx,
			AS, best_i, best_j, out_freqmx);
		asserta(out_freqmx.size() == (AS-1)*(AS-1));
		in_freqmx = out_freqmx;
		ProgressLog("[%3u]  %3u  %3u  H=%.4f\n",
			AS, best_i, best_j, H);

		upd_map(fullcode2reducedcode, best_i, best_j);
		}
	}

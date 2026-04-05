#include "myutils.h"
#include "flat_helpers.h"

static uint get_maxreducedcode(
	const vector<uint> &fullcode2reducedcode)
	{
	const uint fullAS = uint(fullcode2reducedcode.size());
	uint maxreducedcode = 0;
	for (uint full_code = 0; full_code < fullAS ; ++full_code)
		{
		uint reduced_code = fullcode2reducedcode[full_code];
		maxreducedcode = max(reduced_code, maxreducedcode);
		}
	asserta(maxreducedcode > 0);
	return maxreducedcode;
	}

/** After merge, merged pair -> newAS-1; other old reduced codes -> 0..newAS-2. */
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

	//vector<double> out_freqmx;
	vector<uint> fullcode2reducedcode;

	for (uint i = 0; i < in_alpha_size; ++i)
		fullcode2reducedcode.push_back(i);

	validate_map(fullcode2reducedcode);
	for (uint AS = in_alpha_size; AS >= out_alpha_size; --AS)
		{
		uint best_i, best_j;
		double H = find_best_ij(in_freqmx,
			AS, best_i, best_j, fullcode2reducedcode);
		ProgressLog("[%3u]  %3u  %3u  H=%.4f\n",
			AS, best_i, best_j, H);

		upd_map(fullcode2reducedcode, best_i, best_j);
		}
	}

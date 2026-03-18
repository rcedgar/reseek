#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "tabbedlines.h"

void trunc_label(const string &Label,
	string &TruncatedLabel);

void read_fasta_label2idx(
	const string &fafn,
	map<string, uint> &label2idx)
	{
	label2idx.clear();
	SeqDB db_fa;
	db_fa.FromFasta(fafn, false);
	const uint nseqs = db_fa.GetSeqCount();
	for (uint i = 0; i < nseqs; ++i)
		{
		string label;
		trunc_label(db_fa.GetLabel(i), label);
		label2idx[label] = i;
		}
	}

void read_feature_fasta(
	const string &fafn,
	uint alpha_size,
	const map<string, uint> &label2idx,
	vector<vector<uint8_t> > &codeseqs)
	{
	codeseqs.clear();

	const uint8_t *char2letter = (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);

	SeqDB db_fa;
	db_fa.FromFasta(fafn, false);
	const uint nseqs = db_fa.GetSeqCount();
	asserta(nseqs == SIZE(label2idx));
	codeseqs.resize(nseqs);
	for (uint i = 0; i < nseqs; ++i)
		{
		uint L = db_fa.GetSeqLength(i);
		string label;
		trunc_label(db_fa.GetLabel(i), label);
		map<string, uint>::const_iterator iter = label2idx.find(label);
		if (iter == label2idx.end())
			Die("Not found %s in %s", label.c_str(), fafn.c_str());
		uint idx = iter->second;
		asserta(idx < nseqs);
		asserta(codeseqs[idx].size() == 0);
		codeseqs[idx].resize(L);
		const byte *byteseq = db_fa.GetByteSeq(i);
		for (uint pos = 0; pos < L; ++pos)
			{
			uint8_t code = char2letter[byteseq[pos]];
			if (code < alpha_size)
				codeseqs[idx][pos] = code;
			else if (code == 0xff)
				codeseqs[idx][pos] = 0;
			else
				asserta(false);
			}
		}
	}

void read_profiles_from_fastas(
	const vector<string> &fafns,
	const vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles)
	{
	profiles.clear();
	labels.clear();

	const uint nfeat = SIZE(fafns);
	asserta(SIZE(alpha_sizes) == nfeat);
	asserta(nfeat > 0);

	map<string, uint> label2idx;
	read_fasta_label2idx(fafns[0], label2idx);
	const uint nprof = SIZE(label2idx);

	vector<vector<vector<uint8_t> > > codeseqsvec(nfeat);

	read_feature_fasta(fafns[0], alpha_sizes[0], label2idx, codeseqsvec[0]);
	for (uint fi = 1; fi < nfeat; ++fi)
		read_feature_fasta(fafns[fi], alpha_sizes[fi], label2idx, codeseqsvec[fi]);

	profiles.resize(nprof);
	for (auto iter : label2idx)
		{
		void TruncLabel(string &lab);
		string label = iter.first;
		TruncLabel(label);
		uint idx = iter.second;
		labels.push_back(label);
		vector<uint8_t> &profile = profiles[idx];
		uint L = SIZE(codeseqsvec[0][idx]);
		profile.resize(nfeat*L, 0xff);
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			const vector<uint8_t> &codeseq = codeseqsvec[fi][idx];
			for (uint i = 0; i < L; ++i)
				asserta(codeseq[i] != 0xff);
			for (uint pos = 0; pos < L; ++pos)
				profile[fi*L + pos] = codeseq[pos];
			}
		for (uint i = 0; i < nfeat*L; ++i)
			asserta(profile[i] != 0xff);
		}
	}

uint lines2logoddsmx(
	const vector<string> &lines,
	vector<float> &logoddsmx)
	{
	logoddsmx.clear();
	tabbedlines tl(lines);
	uint alpha_size = tl.get_int("logodds");
	asserta(alpha_size != 0);
	logoddsmx.resize(alpha_size*alpha_size);
	tl.get_float_flat_square_mx(alpha_size, logoddsmx.data());
	return alpha_size;
	}

uint read_logoddsmx(
	const string &fn,
	vector<float> &logoddsmx)
	{
	vector<string> lines;
	ReadLinesFromFile(fn, lines);
	return lines2logoddsmx(lines, logoddsmx);
	}

void read_logoddsvec(
	const vector<string> &fns,
	vector<vector<float> > &logoddsmxvec)
	{
	logoddsmxvec.clear();
	for (auto fn : fns)
		{
		vector<float> logoddsmx;
		read_logoddsmx(fn, logoddsmx);
		logoddsmxvec.push_back(logoddsmx);
		}
	}

void log_flat_square_mx(const float *mx, uint n)
	{
	for (uint i = 0; i < n; ++i)
		{
		Log("%2u  |", i);
		for (uint j = 0; j < n; ++j)
			Log(" %7.3g", mx[n*i + j]);
		Log("\n");
		}
	}

void check_profile(
	vector<uint8_t> &profile,
	vector<uint> &alpha_sizes)
	{
	const uint n = SIZE(profile);
	const uint nfeat = SIZE(alpha_sizes);
	asserta(n%nfeat == 0);
	const uint L = n/nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = alpha_sizes[fi];
		for (uint i = 0; i < L; ++i)
			asserta(profile[fi*L + i] < AS);
		}
	}

void check_profiles(
	vector<vector<uint8_t> > &profiles,
	vector<uint> &alpha_sizes)
	{
	const uint nprof = SIZE(profiles);
	for (uint i = 0; i < nprof; ++i)
		check_profile(profiles[i], alpha_sizes);
	}

void read_profiles_and_logoddsmxvec(
	const string &specfn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles,
	vector<vector<float> > &logoddsmxvec)
	{
	vector<string> lines;
	ReadLinesFromFile(specfn, lines);
	uint nfeat = SIZE(lines);

	vector<string> fafns;
	vector<string> logoddsfns;
	vector<float> weights;
	float sumw = 0;
	for (auto line : lines)
		{
		vector<string> flds;
		Split(line, flds, '\t');
		asserta(SIZE(flds) == 5);
		feature_names.push_back(flds[0]);
		fafns.push_back(flds[1]);
		alpha_sizes.push_back(StrToUint(flds[2]));
		logoddsfns.push_back(flds[3]);
		float w = (float) StrToFloat(flds[4]);
		weights.push_back(w);
		sumw += w;
		}
	asserta(SIZE(weights) == nfeat);

	asserta(sumw > 0);
	float sumw2 = 0;
	for (uint i = 0; i < nfeat; ++i)
		{
		float w = weights[i]/sumw;
		sumw2 += w;
		weights[i] = w;
		}
	asserta(sumw2 > 0.99 && sumw2 < 1.01);

	read_logoddsvec(logoddsfns, logoddsmxvec);
	asserta(SIZE(logoddsmxvec) == nfeat);

	read_profiles_from_fastas(fafns, alpha_sizes, labels, profiles);
	}

void cmd_flat_profiles()
	{
	const string &specfn = g_Arg1;
	vector<string> labels;
	vector<vector<uint8_t> > profiles;
	vector<string> feature_names;
	vector<uint> alpha_sizes;
	vector<vector<float> > logoddsmxvec;
	read_profiles_and_logoddsmxvec(
		specfn,
		feature_names,
		alpha_sizes,
		labels,
		profiles,
		logoddsmxvec);

	const uint nfeat = SIZE(feature_names);
	const uint nprof = SIZE(labels);

	asserta(SIZE(alpha_sizes) == nfeat);
	asserta(SIZE(profiles) == nprof);

	check_profiles(profiles, alpha_sizes);

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		Log("\n%s\n", feature_names[fi].c_str());
		log_flat_square_mx(logoddsmxvec[fi].data(), alpha_sizes[fi]);
		}

	// Convert profiles back to FASTA for correctness checking
	if (optset_output2)
		{
		const uint nprof = SIZE(labels);
		asserta(SIZE(profiles) == nprof);
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			uint alpha_size = alpha_sizes[fi];
			const uint8_t *letter2char = (alpha_size == 20 ? g_LetterToCharAmino : g_LetterToCharMu);

			string fn = opt(output2) + feature_names[fi];
			Progress("FASTA %s\n", fn.c_str());
			FILE *ffa = CreateStdioFile(fn);
			for (uint seqidx = 0; seqidx < nprof; ++seqidx)
				{
				const string &label = labels[seqidx];
				const vector<uint8_t> &profile = profiles[seqidx];
				asserta(SIZE(profile)%nfeat == 0);
				const uint L = SIZE(profile)/nfeat;
				string seq;
				for (uint pos = 0; pos < L; ++pos)
					{
					uint8_t code = profile[fi*L + pos];
					seq += letter2char[code];
					}
				SeqToFasta(ffa, label, seq);
				}
			CloseStdioFile(ffa);
			}
		}
	}

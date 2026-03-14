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
	const uint nseq = SIZE(label2idx);

	vector<vector<vector<uint8_t> > > codeseqsvec(nfeat);

	read_feature_fasta(fafns[0], alpha_sizes[0], label2idx, codeseqsvec[0]);
	for (uint fi = 1; fi < nfeat; ++fi)
		read_feature_fasta(fafns[fi], alpha_sizes[fi], label2idx, codeseqsvec[fi]);

	profiles.resize(nseq);
	for (auto iter : label2idx)
		{
		const string &label = iter.first;
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
				profile[nfeat*pos + fi] = codeseq[pos];
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

void cmd_flat_profiles()
	{
	const string &specfn = g_Arg1;
	vector<string> lines;
	ReadLinesFromFile(specfn, lines);
	uint nfeat = SIZE(lines);

	vector<string> fafns;
	vector<string> logoddsfns;
	vector<uint> alpha_sizes;
	for (auto line : lines)
		{
		vector<string> flds;
		Split(line, flds, '\t');
		asserta(SIZE(flds) == 3);
		fafns.push_back(flds[0]);
		alpha_sizes.push_back(StrToUint(flds[1]));
		logoddsfns.push_back(flds[2]);
		}

	vector<string> labels;
	vector<vector<uint8_t> > profiles;
	read_profiles_from_fastas(fafns, alpha_sizes, labels, profiles);

	vector<vector<float> > logoddsmxvec;
	read_logoddsvec(logoddsfns, logoddsmxvec);
	asserta(SIZE(logoddsmxvec) == nfeat);

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		Log("\n%s\n", fafns[fi].c_str());
		log_flat_square_mx(logoddsmxvec[fi].data(), alpha_sizes[fi]);
		}
	}

#include "myutils.h"
#include "flat_profiles.h"
#include "flat_features.h"
#include "flat_helpers.h"
#include "seqdb.h"

void flat_profiles::read_profiles_from_fastas(
	const vector<string> &fafns,
	const unordered_map<string, uint> &label2idx)
	{
	m_label2idx = label2idx;
	const uint nfeat = flat_features::m_nfeat;
	asserta(nfeat);
	asserta(fafns.size() == nfeat);
	const uint32_t *alpha_sizes = flat_features::m_alpha_sizes;

	vector<vector<vector<uint8_t> > > codeseqsvec(nfeat);

	for (uint fi = 0; fi < nfeat; ++fi)
		read_feature_fasta(fafns[fi], alpha_sizes[fi], label2idx, codeseqsvec[fi]);
	const uint nprof = uint(codeseqsvec[0].size());
	for (uint fi = 0; fi < nfeat; ++fi)
		asserta(codeseqsvec[fi].size() == nprof);

	m_profiles.resize(nprof);
	m_labels.resize(nprof);
	for (auto iter : label2idx)
		{
		string label = iter.first;
		uint idx = iter.second;
		m_labels[idx] = label;
		vector<uint8_t> &profile = m_profiles[idx];
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

void flat_profiles::read_profiles_faprof(
	const string &fn,
	vector<string> &feature_names)
	{
	asserta(fn != "");
	asserta(m_labels.empty());
	asserta(m_profiles.empty());
	m_label2idx.clear();
	feature_names.clear();

	m_labels.clear();
	m_profiles.clear();

	SeqDB DB;
	DB.FromFasta(fn);
	const uint ndbseq = DB.GetSeqCount();
	asserta(ndbseq > 0);
	vector<string> flds;
	vector<string> flds2;
	string acc;
	uint nfeat = 0;
	for (uint i = 0; i < ndbseq; ++i)
		{
		const string &label = DB.GetLabel(i);
		Split(label, flds, ':');
		asserta(SIZE(flds) == 2);
		if (i == 0)
			{
			acc = flds[0];
			continue;
			}
		else
			{
			if (flds[0] != acc)
				{
				nfeat = i;
				break;
				}
			}
		}
	asserta(nfeat > 0);
	for (uint i = 0; i < nfeat; ++i)
		{
		const string &label = DB.GetLabel(i);
		Split(label, flds, ':');
		asserta(flds.size() == 2);
		feature_names.push_back(flds[1]);
		}
	ProgressLog("%u features", nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		ProgressLog(" %s", feature_names[fi].c_str());
	ProgressLog("\n");
	asserta(ndbseq%nfeat == 0);
	uint nprof = ndbseq/nfeat;
	m_profiles.resize(nprof);
	for (uint profidx = 0; profidx < nprof; ++profidx)
		{
		uint L = DB.GetSeqLength(nfeat*profidx);
		uint profile_length = nfeat*L;
		vector<uint8_t> &profile = m_profiles[profidx];
		profile.resize(profile_length);
		uint k = 0;
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			const uint seqidx = nfeat*profidx + fi;
			const string &label = DB.GetLabel(seqidx);

			Split(label, flds, ':');
			asserta(flds.size() == 2);
			const string &acc = flds[0];
			if (fi == 0)
				{
				m_labels.push_back(acc);
				m_label2idx[acc] = profidx;
				}
			else
				asserta(acc == m_labels.back());
			const uint8_t *char2letter =
				get_char2letter(feature_names[fi]);
			uint L_fi = DB.GetSeqLength(seqidx);
			asserta(L_fi == L);
			const string &seq = DB.GetSeq(seqidx);
			asserta(SIZE(seq) == L);
			for (uint k = 0; k < L; ++k)
				profile[fi*L + k] = char2letter[seq[k]];
			}
		}
	}

void flat_profiles::profile_to_fasta(FILE *f, uint i) const
	{
	if (f == 0)
		return;
	asserta(i < m_profiles.size());
	asserta(i < m_labels.size());
	const vector<uint8_t> &profile = m_profiles[i];
	const string &label = m_labels[i];
	const uint n = SIZE(profile);
	const uint nfeat = flat_features::m_nfeat;
	asserta(n%nfeat == 0);
	const uint L = n/nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint alpha_size = flat_features::m_alpha_sizes[fi];
		const string &feature_name = flat_features::m_feature_names[fi];
		const uint8_t *letter2char = get_letter2char(alpha_size);

		string seq;
		for (uint pos = 0; pos < L; ++pos)
			{
			uint8_t code = profile[fi*L + pos];
			seq += letter2char[code];
			}
		string label_feat;
		Psa(label_feat, "%s:%s",
			label.c_str(),
			feature_name.c_str());
		SeqToFasta(f, label_feat, seq, L);
		}
	}

void flat_profiles::check_profile(uint i) const
	{
	asserta(i < m_profiles.size());
	const vector<uint8_t> &profile = m_profiles[i];
	const uint n = SIZE(profile);
	const uint nfeat = flat_features::m_nfeat;
	asserta(n%nfeat == 0);
	const uint L = n/nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = flat_features::m_alpha_sizes[fi];
		for (uint i = 0; i < L; ++i)
			asserta(profile[fi*L + i] < AS);
		}
	}

void flat_profiles::check_profiles() const
	{
	const uint nprof = get_nprof();
	for (uint i = 0; i < nprof; ++i)
		check_profile(i);
	}

uint8_t *flat_profiles::get_rev_profile(uint i) const
	{
	asserta(i < m_profiles.size());
	const vector<uint8_t> &profile = m_profiles[i];
	const uint n = SIZE(profile);
	const uint nfeat = flat_features::m_nfeat;
	asserta(n%nfeat == 0);
	const uint L = n/nfeat;
	uint8_t *rev_profile = myalloc(uint8_t, n);
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = flat_features::m_alpha_sizes[fi];
		for (uint i = 0; i < L; ++i)
			{
			uint rev_i = L - i - 1;
			rev_profile[fi*L + rev_i] = profile[fi*L + i];
			}
		}
	return rev_profile;
	}

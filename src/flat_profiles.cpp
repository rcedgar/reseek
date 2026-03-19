#include "myutils.h"
#include "flat_profiles.h"
#include "flat_features.h"
#include "flat_helpers.h"
#include "seqdb.h"

void flat_profiles::read_profiles_faprof(
	const string &fn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes)
	{
	asserta(m_ff == 0); // must create later
	asserta(fn != "");
	asserta(m_labels.empty());
	asserta(m_profiles.empty());
	feature_names.clear();
	alpha_sizes.clear();

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
		string alpha_annot = flds[1];
		Split(alpha_annot, flds, '*');
		asserta(flds.size() == 2);
		feature_names.push_back(flds[0]);
		alpha_sizes.push_back(StrToUint(flds[1]));
		}
	ProgressLog("%u features", nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		ProgressLog(" %s*%u",
			feature_names[fi].c_str(),
			alpha_sizes[fi]);
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
			const string &label = DB.GetLabel(nfeat*profidx + fi);

			Split(label, flds, ':');
			asserta(flds.size() == 2);
			const string &acc = flds[0];
			string name_star_size = flds[1];

			Split(name_star_size, flds2, '*');
			asserta(flds2.size() == 2);
			asserta(flds2[0] == feature_names[fi]);
			asserta(StrToUint(flds2[1]) == alpha_sizes[fi]);
			if (fi == 0)
				m_labels.push_back(acc);
			else
				asserta(acc == m_labels.back());
			uint alpha_size = alpha_sizes[fi];
			const uint8_t *char2letter = get_char2letter(alpha_size);
			uint L_fi = DB.GetSeqLength(nfeat*profidx + fi);
			asserta(L_fi == L);
			const string &seq = DB.GetSeq(nfeat*profidx + fi);
			asserta(SIZE(seq) == L);
			for (uint k = 0; k < L; ++k)
				profile[fi*L + k] = char2letter[seq[k]];
			}
		}
	asserta(m_ff == 0); // must create later
	}

void flat_profiles::check_profile(uint i) const
	{
	asserta(m_ff != 0);
	asserta(i < m_profiles.size());
	const vector<uint8_t> &profile = m_profiles[i];
	const uint n = SIZE(profile);
	const uint nfeat = m_ff->m_nfeat;
	asserta(n%nfeat == 0);
	const uint L = n/nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = m_ff->m_alpha_sizes[fi];
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

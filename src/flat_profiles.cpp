#include "myutils.h"
#include "flat_profiles.h"
#include "flat_alphas.h"
#include "flat_helpers.h"
#include "seqdb.h"

void flat_profiles::read_profiles_from_fastas(
	const vector<string> &fafns,
	const unordered_map<string, uint> &label2idx)
	{
	asserta(m_labels.empty());
	asserta(m_lengths.empty());
	asserta(m_profiles.empty());
	asserta(m_label2idx.empty());

	m_label2idx = label2idx;
	const uint nfeat = flat_alphas::m_nfeat;
	asserta(nfeat);
	asserta(fafns.size() == nfeat);
	const uint32_t *alpha_sizes = flat_alphas::m_alpha_sizes;

	vector<vector<vector<uint8_t> > > codeseqsvec(nfeat);

	for (uint fi = 0; fi < nfeat; ++fi)
		read_feature_fasta(fafns[fi], alpha_sizes[fi], label2idx, codeseqsvec[fi]);
	const uint nprof = uint(codeseqsvec[0].size());
	for (uint fi = 0; fi < nfeat; ++fi)
		asserta(codeseqsvec[fi].size() == nprof);

	m_profiles.resize(nprof);
	m_labels.resize(nprof);
	m_lengths.resize(nprof);
	for (auto iter : label2idx)
		{
		string label = iter.first;
		uint idx = iter.second;
		//vector<uint8_t> &profile = m_profiles[idx];
		uint L = SIZE(codeseqsvec[0][idx]);
		uint8_t *profile = myalloc(uint8_t, nfeat*L);
		memset(profile, 0xff, nfeat*L);
		//profile.resize(nfeat*L, 0xff);
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

		m_labels[idx] = label;
		m_profiles[idx] = profile;
		m_lengths[idx] = L;
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
		uint8_t *profile = myalloc(uint8_t, nfeat*L);
		m_profiles[profidx] = profile;
		//profile.resize(profile_length);
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
				m_labels.push_back(acc.c_str());
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
	const uint8_t *profile = m_profiles[i];
	const string &label = m_labels[i];
	const uint L = m_lengths[i];
	const uint nfeat = flat_alphas::m_nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint alpha_size = flat_alphas::m_alpha_sizes[fi];
		const string &alpha_name = flat_alphas::m_alpha_names[fi];
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
			alpha_name.c_str());
		SeqToFasta(f, label_feat, seq, L);
		}
	}

void flat_profiles::check_profile(uint i) const
	{
	asserta(i < m_profiles.size());
	const uint8_t *profile = m_profiles[i];
	const uint L = m_lengths[i];
	const uint nfeat = flat_alphas::m_nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = flat_alphas::m_alpha_sizes[fi];
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
	const uint8_t *profile = m_profiles[i];
	const uint L = m_lengths[i];
	const uint nfeat = flat_alphas::m_nfeat;
	uint8_t *rev_profile = myalloc(uint8_t, L*nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = flat_alphas::m_alpha_sizes[fi];
		for (uint i = 0; i < L; ++i)
			{
			uint rev_i = L - i - 1;
			rev_profile[fi*L + rev_i] = profile[fi*L + i];
			}
		}
	return rev_profile;
	}

void flat_profiles::from_chains_lookup(const lookup &look,
	const vector<flat_chain_t *> &chains)
	{
	const uint ndom = look.get_ndom();
	const uint nfeat = flat_alphas::m_nfeat;
	asserta(nfeat > 0);

	// Allocate in lookup/domidx order so all downstream code can index by domidx.
	m_labels.clear();
	m_lengths.clear();
	m_profiles.clear();
	m_label2idx.clear();

	m_labels.resize(ndom);
	m_lengths.resize(ndom, 0);
	m_profiles.resize(ndom, 0);

	vector<bool> found(ndom, false);

	for (uint i = 0; i < uint(chains.size()); ++i)
		{
		const flat_chain_t *chain = chains[i];
		if (chain == 0)
			continue;

		string label = chain->m_label;
		trunc_label(label);
		uint domidx = look.get_domidx(label, true);
		if (domidx == UINT_MAX)
			continue; // chain not in lookup: ignore
		m_profiles[domidx] = make_profile(*chain);
		found[domidx] = true;
		m_lengths[domidx] = chain->get_length();
		m_labels[domidx] = label;
		}

	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		if (!found[domidx])
			Die("Missing chain >%s",
				look.get_dom(domidx).c_str());
		}

	check_profiles();
	}

void flat_profiles::from_chains(const vector<flat_chain_t *> &chains)
	{
	asserta(m_labels.empty());
	asserta(m_profiles.empty());
	asserta(m_label2idx.empty());

	const uint nchain = uint(chains.size());
	const uint nfeat = flat_alphas::m_nfeat;
	asserta(nfeat);
	const uint32_t *alpha_sizes = flat_alphas::m_alpha_sizes;

	m_profiles.resize(nchain);
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		const flat_chain_t &chain = *chains[chainidx];
		m_labels.push_back(chain.m_label);
		const uint L = chain.get_length();
		if (L == 0) continue;
		m_profiles[chainidx] = make_profile(chain);
		}
	}

uint8_t *flat_profiles::make_profile(const flat_chain_t &chain) const
	{
	const uint L = chain.get_length();
	asserta(L > 0);

	const uint nfeat = flat_alphas::m_nfeat;
	asserta(nfeat > 0);

	uint8_t *profile = myalloc(uint8_t, nfeat*L);
	memset(profile, 0xff, nfeat*L);

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const FAN fan = flat_alphas::m_fans[fi];
		const uint alpha_size = flat_alphas::m_alpha_sizes[fi];

		uint8_t *codeseq = profile + fi*L;
		if (is_quantized(fan))
			chaq::slow_get_codeseq_binned(&chain, fan, alpha_size, codeseq);
		else
			{
			const uint8_t undef_code = chaq::get_undef_code(fan, alpha_size);
			chaq::slow_get_codeseq_discrete(
				&chain, fan, alpha_size, undef_code, codeseq);
			}

#if DEBUG
		for (uint pos = 0; pos < L; ++pos)
			assert(codeseq[pos] < alpha_size);
#endif
		}

	return profile;
	}

void flat_profiles::write_nu_hexfasta(const string &fn) const
	{
	if (fn == "") return;
	FILE *f = CreateStdioFile(fn);
	const uint nprof = uint(m_nu_codeseqs.size());
	for (uint i = 0; i < nprof; ++i)
		{
		const uint8_t *codeseq = m_nu_codeseqs[i];
		const uint L = get_length(i);
		string hexseq;
		hexseq.reserve(L);
		for (uint pos = 0; pos < L; ++pos)
			Psa(hexseq, "%02x", codeseq[pos]);
		const string &label = m_labels[i];
		SeqToFasta(f, label, hexseq);
		}
	CloseStdioFile(f);
	}

void flat_profiles::set_nu_codeseqs(const string &hexfastafn)
	{
	uint fi_aa20 = flat_alphas::get_fi(FAN_aa, 20);
	uint fi_pm2 = flat_alphas::get_fi(FAN_pm, 2);
	uint fi_sec32 = flat_alphas::get_fi(FAN_sec, 32);
	uint nprof = get_nprof();
	m_nu_codeseqs.resize(nprof, 0);
	m_nu_codeseqs_rev.resize(nprof, 0);
	for (uint i = 0; i < nprof; ++i) 
		set_nu_codeseq(fi_aa20, fi_pm2, fi_sec32, i);
	write_nu_hexfasta(hexfastafn);
	}

void flat_profiles::set_nu_codeseq(
	uint fi_aa20, uint fi_pm2, uint fi_sec32, uint idx)
	{
	assert(idx < m_profiles.size());
	const uint8_t *profile = m_profiles[idx];
	uint L = get_length(idx);
	uint8_t *codeseq = myalloc(uint8_t, L);
	uint8_t *codeseq_rev = myalloc(uint8_t, L);
	const uint8_t *prof_aa20 = profile + size_t(fi_aa20)*L;
	const uint8_t *prof_pm2 = profile + size_t(fi_pm2)*L;
	const uint8_t *prof_sec32 = profile + size_t(fi_sec32)*L;

	for (uint pos = 0; pos < L; ++pos)
		{
		uint8_t code_aa20 = prof_aa20[pos];
		uint8_t code_pm2 = prof_pm2[pos];
		uint8_t code_sec32 = prof_sec32[pos];
		assert(code_aa20 < 20);
		assert(code_pm2 < 2);
		assert(code_sec32 < 32);
		uint8_t code_aa4 = chaq::m_aacode2aa4code[code_aa20];
#if DEBUG
		uint32 code_nu = code_aa4 + code_pm2*4 + code_sec32*4*2;
		assert(code_nu < 256);
		codeseq[pos] = int8_t(code_nu);
		codeseq_rev[L-pos-1] = int8_t(code_nu);
#else
		uint8_t code_nu = code_aa4 + code_pm2*4 + code_sec32*4*2;
		codeseq[pos] = code_nu;
		codeseq_rev[L-pos-1] = int8_t(code_nu);
#endif
		}
	m_nu_codeseqs[idx] = codeseq;
	m_nu_codeseqs_rev[idx] = codeseq_rev;
	}

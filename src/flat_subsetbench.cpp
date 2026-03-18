#include "myutils.h"
#include "flat_subsetbench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include <unordered_map>
#include <unordered_set>

void read_profiles_and_logoddsmxvec(
	const string &specfn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles,
	vector<vector<float> > &logoddsmxvec);

uint32_t get_flat_pssm_feature_block_offsets(
	const uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	uint32_t * __restrict feature_block_offsets);

void fill_flat_pssm(
	const uint8_t * __restrict profQ,
	uint32_t LQ,
	uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	const uint32_t * __restrict feature_block_offsets,
	const float *const * __restrict weighted_logoddsmxvec,
	float * __restrict pssm);

float sw_flat_pssm(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profA, uint LA,
	const float *__restrict pssm, uint LB,
	const uint32_t * __restrict feature_block_offsets,
	uint nfeat,
	float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	string &Path);

static const uint MAGIC1 = 0xd06e1;
static const uint MAGIC2 = 0xd06e2;
static const uint MAGIC3 = 0xd06e3;
static const uint MAGIC4 = 0xd06e4;

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void flat_subsetbench::validate_mappings() const
	{
	//Bug("start validate_mappings()");
	asserta(!m_DopeDomIdxs.empty());
	const uint ndom = SIZE(m_Doms);
	for (auto DomIdx : m_DopeDomIdxs)
		{
		asserta(DomIdx < ndom);
		string Dom = m_Doms[DomIdx];
		TruncLabel(Dom);
		map<string, uint>::const_iterator iter =
			m_DomToIdx.find(Dom);
		asserta(iter != m_DomToIdx.end());
		asserta(iter->second == DomIdx);
		}

	for (uint i = 0; i < m_DopeSize; ++i)
		{
		uint DomIdxQ = m_DomIdxQs[i];
		uint DomIdxT = m_DomIdxTs[i];

		asserta(DomIdxQ < ndom);
		asserta(DomIdxT < ndom);
		}

	uint pair_count = 0;
	const uint nuq = SIZE(m_DomIdxQs_sorted);
	asserta(SIZE(m_TargetDomIdxs) == nuq);
	asserta(SIZE(m_TargetPairIdxs) == nuq);
	set<pair<uint16_t, uint16_t> > done_pairs;
	for (uint i = 0; i < nuq; ++i)
		{
		uint16_t DomIdxQ = m_DomIdxQs_sorted[i];
		asserta(DomIdxQ < ndom);
		const vector<uint16_t> &DomIdxs = m_TargetDomIdxs[i];
		const vector<uint16_t> &PairIdxs = m_TargetPairIdxs[i];
		const uint nt = SIZE(DomIdxs);
		asserta(SIZE(PairIdxs) == nt);
		pair_count += nt;
		for (uint ti = 0; ti < nt; ++ti)
			{
			uint16_t DomIdxT = DomIdxs[ti];
			asserta(DomIdxT < ndom);
			pair<uint16_t, uint16_t> pairQT(DomIdxQ, DomIdxT);
			pair<uint16_t, uint16_t> pairTQ(DomIdxT, DomIdxQ);
			asserta(done_pairs.find(pairQT) == done_pairs.end());
			asserta(done_pairs.find(pairTQ) == done_pairs.end());
			done_pairs.insert(pairQT);
			done_pairs.insert(pairTQ);
			}
		}
	asserta(pair_count == m_DopeSize);
	const uint nprof = SIZE(m_profiles);
	const uint nfeat = SIZE(m_alpha_sizes);
	uint missing = 0;
	//Bug("before test loop");
	for (auto DomIdx : m_DopeDomIdxs)
		{
		string Dom = m_Doms[DomIdx];
		TruncLabel(Dom);
		uint profile_idx = m_DomIdx_to_profile_idx[DomIdx];
		if (profile_idx == UINT16_MAX)
			{
			++missing;
			Log("missing >%s\n", Dom.c_str());
			continue;
			}
		asserta(profile_idx < nprof);
		string prof_label = m_profile_labels[profile_idx];
		TruncLabel(prof_label);
		asserta(prof_label == Dom);
		uint profile_length = SIZE(m_profiles[profile_idx]);
		asserta(profile_length%nfeat == 0);
		}
	if (missing != 0) Die("%u missing", missing);

	ProgressLog("validate_mappings() PASSED\n");
	}

void flat_subsetbench::LoadAlphas(const string &SpecFN)
	{
	read_profiles_and_logoddsmxvec(
		SpecFN,
		m_AlphaNames,
		m_alpha_sizes,
		m_profile_labels,
		m_profiles,
		m_raw_logoddsmxvec);

	uint nfeat = SIZE(m_AlphaNames);
	asserta(SIZE(m_alpha_sizes) == nfeat);
	m_weighted_logoddsmxvec.clear();
	m_weighted_logoddsmxvec.resize(nfeat);

// Initialize to uniform weights
	const float w = 1.0f/nfeat;;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = m_alpha_sizes[fi];
		m_weighted_logoddsmxvec[fi].reserve(AS);
		asserta(SIZE(m_raw_logoddsmxvec[fi]) == AS*AS);
		for (uint k = 0; k < AS*AS; ++k)
			m_weighted_logoddsmxvec[fi].push_back(
				w*m_raw_logoddsmxvec[fi][k]);
		}

	const uint n = SIZE(m_profile_labels);
	asserta(SIZE(m_profiles) == n);

	set<uint16_t> domidxset;
	set<string> domlabelset;
	for (uint i = 0; i < m_DopeSize; ++i)
		{
		uint domidxQ = m_DomIdxQs[i];
		uint domidxT = m_DomIdxTs[i];

		domidxset.insert(domidxQ);
		domidxset.insert(domidxT);

		domlabelset.insert(m_Doms[domidxQ]);
		domlabelset.insert(m_Doms[domidxT]);
		}

	unordered_map<string, uint> profile_label_to_idx;
	m_MaxL = 0;
	for (uint profile_idx = 0; profile_idx < n; ++profile_idx)
		{
		const string &label = m_profile_labels[profile_idx];
		profile_label_to_idx[label] = profile_idx;

		uint Ln = SIZE(m_profiles[profile_idx]);
		asserta(Ln%nfeat == 0);
		uint L = Ln/nfeat;
		m_MaxL = max(L, m_MaxL);
		}

	m_DomIdx_to_profile_idx.clear();
	m_DomIdx_to_profile_idx.resize(UINT16_MAX, UINT16_MAX);
	set<string> foundset;
	for (auto DomIdx : domidxset)
		{
		string label = m_Doms[DomIdx];
		TruncLabel(label);
		unordered_map<string, uint>::const_iterator iter =
			profile_label_to_idx.find(label);
		if (iter != profile_label_to_idx.end())
			{
			foundset.insert(label);;
			m_DomIdx_to_profile_idx[DomIdx] = iter->second;
			}
		}
	uint nfound = SIZE(foundset);
	uint ndom = SIZE(domlabelset);
	asserta(ndom >= nfound);
	if (nfound < ndom)
		{
		for (auto label : domlabelset)
			{
			if (foundset.find(label) == foundset.end())
				Log("missing >%s\n", label.c_str());
			}
		Die("%u found, %u missing profiles", nfound, ndom - nfound);
		}
	}

void flat_subsetbench::AddDom(const string &Dom, const string &ScopId)
	{
	vector<string> Fields;
	Split(ScopId, Fields, '.');
	asserta(SIZE(Fields) == 4);
	const string SF = Fields[0] + "." + Fields[1] + "." + Fields[2];
	uint SFIdx = UINT_MAX;
	if (m_SFToIdx.find(SF) == m_SFToIdx.end())
		{
		SFIdx = SIZE(m_SFs);
		m_SFs.push_back(SF);
		m_SFToIdx[SF] = SFIdx;
		}
	else
		SFIdx = m_SFToIdx[SF];
	m_SFIdxs.push_back(SFIdx);
	asserta(SFIdx < SIZE(m_SFIdxToSize));
	m_SFIdxToSize[SFIdx] += 1;

	uint DomIdx = UINT_MAX;
	if (m_DomToIdx.find(Dom) != m_DomToIdx.end())
		{
		Die("Duplicate dom >%s", Dom.c_str());
		DomIdx = m_DomToIdx[Dom];
		}
	DomIdx = SIZE(m_Doms);
	m_Doms.push_back(Dom + "/" + SF);
	m_DomToIdx[Dom] = DomIdx;
	}

void flat_subsetbench::ReadLookup(const string &FileName)
	{
	m_SFIdxToSize.clear();
	m_SFIdxToSize.resize(2000, 0);

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	vector<string> Fields2;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		string Dom = Fields[0];
		TruncLabel(Dom);
		const string &ScopId = Fields[1];
		AddDom(Dom, ScopId);
		}
	CloseStdioFile(f);

	m_NT = 0;
	const uint SFCount = SIZE(m_SFs);
	for (uint SFIdx = 0; SFIdx < SFCount; ++SFIdx)
		{
		uint Size = m_SFIdxToSize[SFIdx];
		asserta(Size > 0);
		m_NT += (Size*(Size - 1))/2;
		}
	m_NT *= 2;
	}

void flat_subsetbench::AllocDope(uint DopeSize)
	{
	asserta(m_DopeSize == 0);
	m_DomIdxQs = myalloc(uint16_t, DopeSize);
	m_DomIdxTs = myalloc(uint16_t, DopeSize);
	m_TPs = myalloc(bool, DopeSize);
	m_DopeSize = DopeSize;
	}

void flat_subsetbench::AllocHits()
	{
	asserta(m_Scores == 0);
	asserta(m_ScoreOrder == 0);
	m_Scores = myalloc(float, m_DopeSize);
	m_ScoreOrder = myalloc(uint, m_DopeSize);
	}

uint16_t flat_subsetbench::GetDomIdx(const string &Label, bool ErrOk) const
	{
	string Dom;
	SCOP40Bench::GetDomFromLabel(Label, Dom);
	map<string, uint>::const_iterator iter = m_DomToIdx.find(Dom);
	if (iter == m_DomToIdx.end())
		{
		if (ErrOk)
			return UINT16_MAX;
		Die("flat_subsetbench::GetDomIdx(%s)", Label.c_str());
		}
	uint Idx = iter->second;
	return Idx;
	}

uint flat_subsetbench::GetSFIdx(uint DomIdx) const
	{
	asserta(DomIdx < SIZE(m_SFIdxs));
	return m_SFIdxs[DomIdx];
	}

// 1=Q, 2=T, 3=Evalue
// Keep lower triangle only
void flat_subsetbench::MakeDopeFromHits(const string &FN)
	{
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	vector<string> LabelQs;
	vector<string> LabelTs;
	uint HitCount = 0;
	uint HighEvalue = 0;
	uint GtCount = 0;
	uint SelfCount = 0;
	uint NotFound = 0;
	while (ReadLineStdioFile(f, Line))
		{
		if (++HitCount%100000 == 0)
			Progress("Hits %u\r", HitCount);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) >= 3);
		const string &LabelQ = Fields[0];
		const string &LabelT = Fields[1];
		if (LabelQ == LabelT)
			{
			++SelfCount;
			continue;
			}
		string DomQ, DomT;
		SCOP40Bench::GetDomFromLabel(LabelQ, DomQ);
		SCOP40Bench::GetDomFromLabel(LabelT, DomT);

		uint DomIdxQ = GetDomIdx(DomQ, true);
		uint DomIdxT = GetDomIdx(DomT, true);
		if (DomIdxQ == UINT_MAX || DomIdxT == UINT_MAX)
			{
			++NotFound;
			continue;
			}

		double Evalue = StrToFloat(Fields[2]);
		if (Evalue >= 10)
			{
			++HighEvalue;
			continue;
			}
		if (LabelQ > LabelT)
			{
			++GtCount;
			continue;
			}

		LabelQs.push_back(LabelQ);
		LabelTs.push_back(LabelT);
		}
	uint DopeSize = SIZE(LabelQs);
	asserta(SIZE(LabelTs) == DopeSize);

	ProgressLog("%10u  Total hits\n", HitCount);
	ProgressLog("%10u  High E-value\n", HighEvalue);
	ProgressLog("%10u  Not found\n", NotFound);
	ProgressLog("%10u  Self-hits\n", SelfCount);
	ProgressLog("%10u  Other triangle\n", GtCount);
	ProgressLog("%10u  Hits saved to dope\n", DopeSize);

	AllocDope(DopeSize);

	uint NT = 0;
	uint NF = 0;
	for (uint Idx = 0; Idx < m_DopeSize; ++Idx)
		{
		const string &LabelQ = LabelQs[Idx];
		const string &LabelT = LabelTs[Idx];

		string DomQ, DomT;
		SCOP40Bench::GetDomFromLabel(LabelQ, DomQ);
		SCOP40Bench::GetDomFromLabel(LabelT, DomT);

		uint DomIdxQ = GetDomIdx(DomQ, false);
		uint DomIdxT = GetDomIdx(DomT, false);

		m_DopeDomIdxs.insert(DomIdxQ);
		m_DopeDomIdxs.insert(DomIdxT);

		m_DomIdxQs[Idx] = DomIdxQ;
		m_DomIdxTs[Idx] = DomIdxT;

		uint SFIdxQ = GetSFIdx(DomIdxQ);
		uint SFIdxT = GetSFIdx(DomIdxT);

		if (SFIdxQ == SFIdxT)
			{
			++NT;
			m_TPs[Idx] = true;
			}
		else
			{
			++NF;
			m_TPs[Idx] = false;
			}
		}

	CloseStdioFile(f);

	ProgressLog("%u / %u doms in dope, %u TPs, %u FPs\n",
		SIZE(m_DopeDomIdxs), SIZE(m_Doms), NT, NF);
	}

void flat_subsetbench::WriteDope(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);

	WriteStdioFile(f, &MAGIC1, sizeof(MAGIC1));
	WriteStdioFile(f, &m_DopeSize, sizeof(m_DopeSize));
	WriteStdioFile(f, (void *) m_DomIdxQs, m_DopeSize*sizeof(m_DomIdxQs[0]));
	WriteStdioFile(f, (void *) m_DomIdxTs, m_DopeSize*sizeof(m_DomIdxTs[0]));
	WriteStdioFile(f, (void *) m_TPs, m_DopeSize*sizeof(m_TPs[0]));
	WriteStdioFile(f, &MAGIC2, sizeof(MAGIC2));

	CloseStdioFile(f);
	}

void flat_subsetbench::ReadDope(const string &FN)
	{
	if (FN == "")
		return;
	m_DopeFN = FN;
	FILE *f = OpenStdioFile(FN);

	uint Word;
	ReadStdioFile(f, &Word, sizeof(Word));
	asserta(Word == MAGIC1);

	uint DopeSize;
	ReadStdioFile(f, &DopeSize, sizeof(DopeSize));
	AllocDope(DopeSize);

	ReadStdioFile(f, (void *) m_DomIdxQs, m_DopeSize*sizeof(m_DomIdxQs[0]));
	ReadStdioFile(f, (void *) m_DomIdxTs, m_DopeSize*sizeof(m_DomIdxTs[0]));
	ReadStdioFile(f, (void *) m_TPs, m_DopeSize*sizeof(m_TPs[0]));

	ReadStdioFile(f, &Word, sizeof(Word));
	asserta(Word == MAGIC2);

	CloseStdioFile(f);

	SortDope();
	}

void flat_subsetbench::LoadStats() const
	{
	uint nq = SIZE(m_DomIdxQs_sorted);
	ProgressLog("\n");
	ProgressLog("Dope file         %s\n", m_DopeFN.c_str());
	ProgressLog("Pairs to align    %u\n", m_DopeSize);
	ProgressLog("Queries (cached)  %u\n", nq);
	ProgressLog("Targets/query     %.1f\n", float(m_DopeSize)/nq);
	ProgressLog("Max length        %u\n", m_MaxL);
	}

void flat_subsetbench::SortDope()
	{
	m_DomIdxQs_sorted.clear();
	m_TargetDomIdxs.clear();
	m_TargetPairIdxs.clear();
	m_DopeDomIdxs.clear();

	vector<uint16_t> tmp_DomIdxQs;
	std::unordered_map<uint16_t, uint32_t> count;
	std::unordered_map<uint16_t, std::unordered_set<uint16_t>> nbrs;

	for (size_t k = 0; k < m_DopeSize; ++k)
		{
		const uint16_t i = m_DomIdxQs[k];
		const uint16_t j = m_DomIdxTs[k];
		asserta(i != j);

		++count[i];
		++count[j];

		nbrs[i].insert(j);
		nbrs[j].insert(i);

		m_DopeDomIdxs.insert(i);
		m_DopeDomIdxs.insert(j);
		}

	// Collect all distinct values.
	tmp_DomIdxQs.reserve(count.size());
	for (const auto& kv : count)
		tmp_DomIdxQs.push_back(kv.first);

	// Sort by decreasing abundance, then increasing value for tie-break.
	std::sort(tmp_DomIdxQs.begin(), tmp_DomIdxQs.end(),
		[&](uint16_t a, uint16_t b)
		{
		if (count[a] != count[b])
			return count[a] > count[b];
		return a < b;
		});

	size_t n = tmp_DomIdxQs.size();
	m_TargetDomIdxs.reserve(n);
	m_TargetPairIdxs.reserve(n);
	uint pair_idx = 0;
	set<pair<uint16_t, uint16_t> > done_pairs;
	for (uint16_t DomIdxQ : tmp_DomIdxQs)
		{
		vector<uint16_t> not_done_targets;
		vector<uint16_t> pair_idxs;
		const unordered_set<uint16_t> &nbrset = nbrs[DomIdxQ];
		for (uint16_t DomIdxT : nbrset)
			{
			pair<uint16_t, uint16_t> pairQT(DomIdxQ, DomIdxT);
			pair<uint16_t, uint16_t> pairTQ(DomIdxT, DomIdxQ);
			if (done_pairs.find(pairQT) != done_pairs.end()) continue;
			if (done_pairs.find(pairTQ) != done_pairs.end()) continue;
			done_pairs.insert(pairQT);
			done_pairs.insert(pairTQ);
			not_done_targets.push_back(DomIdxT);
			pair_idxs.push_back(pair_idx++);
			}
		if (!not_done_targets.empty())
			{
			m_DomIdxQs_sorted.push_back(DomIdxQ);
			m_TargetDomIdxs.push_back(not_done_targets);
			m_TargetPairIdxs.push_back(pair_idxs);
			}
		}
	asserta(pair_idx == m_DopeSize);
	}

void flat_subsetbench::ThreadBody(uint ThreadIdx)
	{
	const uint NQ = SIZE(m_DomIdxQs_sorted);
	const uint nfeat = SIZE(m_AlphaNames);
	asserta(SIZE(m_weighted_logoddsmxvec) == nfeat);
	float **weighted_logoddsmxvec = myalloc(float *, nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = m_alpha_sizes[fi];
		asserta(AS >= 2 && AS < 256);
		weighted_logoddsmxvec[fi] = m_weighted_logoddsmxvec[fi].data();
		}
	uint32_t *feature_block_offsets = myalloc(uint32_t, nfeat);
	const uint32_t sum_alpha_sizes =
		get_flat_pssm_feature_block_offsets(nfeat,
			m_alpha_sizes.data(), feature_block_offsets);

	float *pssmQ = myalloc(float, m_MaxL*sum_alpha_sizes);

	float *scratch_rows = myalloc(float, 2*m_MaxL + 2);
	const float **scratch_pssms = myalloc(const float *, nfeat);
	uint8_t *TB = myalloc(uint8_t, m_MaxL*m_MaxL);
	uint Loi, Loj, Leni, Lenj;
	string Path;

	for (;;)
		{
		uint qidx = m_NextQueryIdx++;
		if (qidx >= NQ)
			return;
		/////////////////////////////////////////////////////////
		// Cache query
		/////////////////////////////////////////////////////////
		uint DomIdxQ = m_DomIdxQs_sorted[qidx];
		uint prof_idxQ = m_DomIdx_to_profile_idx[DomIdxQ];
		const vector<uint8_t> &profvecQ = m_profiles[prof_idxQ];
		uint profile_length = SIZE(profvecQ);
		asserta(profile_length%nfeat == 0);
		uint LQ = profile_length/nfeat;
		assert(LQ <= m_MaxL);
		if (SIZE(profvecQ) != LQ*nfeat)
			{
			ProgressLog("qidx       %u\n", qidx);
			ProgressLog("nfeat      %u\n", nfeat);
			ProgressLog("DomIdxQ    %u\n", DomIdxQ);
			ProgressLog("prof_idxq  %u\n", prof_idxQ);
			ProgressLog("dom        %s\n", m_Doms[DomIdxQ].c_str());
			ProgressLog("prof       %s\n", m_profile_labels[prof_idxQ].c_str());
			ProgressLog("LQ         %u\n", LQ);
			ProgressLog("proflen    %u\n", SIZE(profvecQ));
			ProgressLog("LQ*nfeat   %u\n", LQ*nfeat);
			Die("SIZE(profvecQ) != LQ*nfeat");
			}
		const uint8_t *profQ = profvecQ.data();
		fill_flat_pssm(profQ, LQ, nfeat, m_alpha_sizes.data(),
			feature_block_offsets, weighted_logoddsmxvec, pssmQ);
		/////////////////////////////////////////////////////////

		const vector<uint16_t> &DomIdxTs = m_TargetDomIdxs[qidx];
		const vector<uint16_t> &PairIdxs = m_TargetPairIdxs[qidx];
		const uint target_count = SIZE(DomIdxTs);
		for (uint ti = 0; ti < target_count; ++ti)
			{
			uint DomIdxT = DomIdxTs[ti];
			uint PairIdx = PairIdxs[ti];

			uint prof_idxT = m_DomIdx_to_profile_idx[DomIdxT];
			const vector<uint8_t> &profvecT = m_profiles[prof_idxT];
			uint profile_lengthT = SIZE(profvecT);
			assert(profile_lengthT%nfeat == 0);
			uint LT = profile_lengthT/nfeat;
			assert(LT <= m_MaxL);
			const uint8_t *profT = profvecT.data();

			float Score = sw_flat_pssm(
				scratch_rows, TB, scratch_pssms,
				profT, LT,
				pssmQ, LQ,
				feature_block_offsets, nfeat,
				-m_Open, -m_Ext,
				Loi, Loj, Leni, Lenj, Path);

			asserta(!isnan(Score));
			asserta(!isinf(Score));
			m_Scores[PairIdx] = Score;
			}
		}
	}

void flat_subsetbench::Search()
	{
	m_ThreadCount = GetRequestedThreadCount();
	m_NextQueryIdx = 0;
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBody, this, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}

void flat_subsetbench::StaticThreadBody(flat_subsetbench *SB, uint ThreadIdx)
	{
	SB->ThreadBody(ThreadIdx);
	}

void flat_subsetbench::SetScoreOrder()
	{
	if (m_ScoreOrder == 0)
		m_ScoreOrder = myalloc(uint, m_DopeSize);
	QuickSortOrderDesc(m_Scores, m_DopeSize, m_ScoreOrder);
	}

void flat_subsetbench::Bench(const string &Msg)
	{
	SetScoreOrder();
	asserta(m_ScoreOrder != 0);
	asserta(m_NT > 0);
	uint K = m_DopeSize;
	uint nt = 0;
	uint nf = 0;
	float LastScore = FLT_MAX;
	float SEPQ0_1 = FLT_MAX;
	float SEPQ1 = FLT_MAX;
	float SEPQ10 = FLT_MAX;
	const uint DomCount = GetDomCount();
	for (uint k = 0; k < K; ++k)
		{
		uint HitIdx = m_ScoreOrder[k];
		float Score = m_Scores[HitIdx];
		if (Score != LastScore)
			{
			if (Score >= LastScore)
				Die("k=%u HitIdx=%u Score=%.3g LastScore=%3g\n",
					k, HitIdx, Score, LastScore);
			float EPQ = 2*float(nf)/DomCount;
			float Sens = 2*float(nt)/m_NT;
			if (SEPQ0_1 == FLT_MAX && EPQ >= 0.1) SEPQ0_1 = Sens;
			if (SEPQ1 == FLT_MAX   && EPQ >= 1)   SEPQ1   = Sens;
			if (SEPQ10 == FLT_MAX  && EPQ >= 10)  SEPQ10  = Sens;
			LastScore = Score;
			}
		if (m_TPs[HitIdx])
			++nt;
		else
			++nf;
		}
	float EPQ = 2*float(nf)/DomCount;
	float Sens = 2*float(nt)/m_NT;
	if (SEPQ0_1 == FLT_MAX) SEPQ0_1 = Sens;
	if (SEPQ1 == FLT_MAX)   SEPQ1   = Sens;
	if (SEPQ10 == FLT_MAX)  SEPQ10  = Sens;
	m_Sum3 = SEPQ0_1*2 + SEPQ1*3/2 + SEPQ10;

	if (Msg != "")
		ProgressLog("%s ", Msg.c_str());
	ProgressLog("SEPQ0.1=%.3f", SEPQ0_1);
	ProgressLog(" SEPQ1=%.3f", SEPQ1);
	ProgressLog(" SEPQ10=%.3f", SEPQ10);
	ProgressLog(" Sum3=%.3f", m_Sum3);
	ProgressLog("\n");
	}

void flat_subsetbench::WriteHits(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);
	for (uint k = 0; k < m_DopeSize; ++k)
		{
		uint HitIdx = m_ScoreOrder[k];
		float Score = m_Scores[HitIdx];
		uint DomIdxQ = m_DomIdxQs[HitIdx];
		uint DomIdxT = m_DomIdxTs[HitIdx];
		const string &DomQ = m_Doms[DomIdxQ];
		const string &DomT = m_Doms[DomIdxT];
		asserta(DomQ != DomT);
		fprintf(f, "%s", DomQ.c_str());
		fprintf(f, "\t%s", DomT.c_str());
		fprintf(f, "\t%.3g", Score);
		fprintf(f, "\n");

		fprintf(f, "%s", DomT.c_str());
		fprintf(f, "\t%s", DomQ.c_str());
		fprintf(f, "\t%.3g", Score);
		fprintf(f, "\n");
		}

	CloseStdioFile(f);
	}

float flat_subsetbench_AF_SWFast(
	XDPMem &Mem,
	uint LQ, uint LT,
	float Open, float Ext,
	float * const * SWMx)
	{
	float SWFast(XDPMem &Mem, const float * const *SMxData, uint LA, uint LB,
	  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	  string &Path);

	uint Loi, Loj, Leni, Lenj;
	string Path;
	float Score = SWFast(Mem, SWMx, LQ, LT, -Open, -Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

void flat_subsetbench::SetScalarParams(
	const vector<string> &Names,
	const vector<float> &Values)
	{
	m_Open = FLT_MAX;
	m_Ext = FLT_MAX;
	for (uint i = 0; i < SIZE(Names); ++i)
		{
		const string &Name = Names[i];
		float Value = Values[i];
		if (Name == "open")
			m_Open = Value;
		else if (Name == "ext")
			m_Ext = Value;
		else if (Name == "gap2")
			{
			m_Open = Value;
			m_Ext = Value/10;
			}
		else
			Die("flat_subsetbench::SetScalarParams() Name=%s", Name.c_str());
		}
	asserta(m_Open != FLT_MAX && m_Ext != FLT_MAX);
	}

void flat_subsetbench::ClassifyParams(
	const vector<string> &Names,
	const vector<float> &Values,
	vector<string> &AlphaNames,
	vector<float> &Weights,
	vector<string> &ScalarNames,
	vector<float> &ScalarValues)
	{
	for (uint i = 0; i < SIZE(Names); ++i)
		{
		const string &Name = Names[i];
		float Value = Values[i];
		if (Name == "open" || Name == "ext" || Name == "gap2")
			{
			ScalarNames.push_back(Name);
			ScalarValues.push_back(Value);
			}
		else
			{
			AlphaNames.push_back(Name);
			Weights.push_back(Value);
			}
		}
	}

void flat_subsetbench::ApplyWeightsToLogOdds(
	const unordered_map<string, float> &NameToWeight)
	{
	uint nalpha = SIZE(m_AlphaNames);
	m_Weights.clear();
	asserta(SIZE(NameToWeight) == nalpha);
	unordered_map<string, uint> NameToIdx;
	for (uint idx = 0; idx < nalpha; ++idx)
		NameToIdx[m_AlphaNames[idx]] = idx;

	for (unordered_map<string, float>::const_iterator iter = NameToWeight.begin();
		iter != NameToWeight.end(); ++iter)
		{
		const string &Name = iter->first;
		float Weight = iter->second;
		unordered_map<string, uint>::const_iterator iter2 = NameToIdx.find(Name);
		asserta(iter2 != NameToIdx.end());
		uint idx = iter2->second;
		m_Weights[idx] = Weight;

		uint AS = m_alpha_sizes[idx];
		for (uint code = 0; code < AS; ++code)
			m_weighted_logoddsmxvec[idx][code] =
				m_raw_logoddsmxvec[idx][code]*Weight;
		}
	}

void flat_subsetbench::UpdateParamsFromVarStr(const string &VarStr)
	{
	vector<string> Names;
	vector<float> Values;
	ParseVarStr(VarStr, Names, Values);

	vector<string> AlphaNames;
	vector<float> Weights;
	vector<string> ScalarNames;
	vector<float> ScalarValues;
	flat_subsetbench::ClassifyParams(
		Names, Values, AlphaNames, Weights, ScalarNames, ScalarValues);

	SetScalarParams(ScalarNames, ScalarValues);

	uint n = SIZE(AlphaNames);
	asserta(SIZE(Weights) == n);
	unordered_map<string, float> NameToWeight;
	for (uint i = 0; i < n; ++i)
		NameToWeight[AlphaNames[i]] = Weights[i];
	ApplyWeightsToLogOdds(NameToWeight);
	}

void flat_subsetbench::InitFB()
	{
	m_FB.m_Labels = m_Doms;
	m_FB.SetLookupFromLabels();
	m_FB.Alloc();
	}

void cmd_flat_subset_bench()
	{
	asserta(optset_spec);
	asserta(optset_lookup);
	asserta(optset_varstr);

	const string &DopeFN = g_Arg1;
	const string &LookupFN = opt(lookup);
	const string &SpecFN = opt(spec);
	const string &VarStr = opt(varstr);

	vector<string> Names;
	vector<float> Values;
	ParseVarStr(VarStr, Names, Values);

	vector<string> AlphaNames;
	vector<float> Weights;
	vector<string> ScalarNames;
	vector<float> ScalarValues;
	flat_subsetbench::ClassifyParams(
		Names, Values, AlphaNames, Weights, ScalarNames, ScalarValues);

	flat_subsetbench SB;
	SB.ReadLookup(LookupFN);
	SB.ReadDope(DopeFN);
	SB.LoadAlphas(SpecFN);
	SB.LoadStats();
	SB.validate_mappings();
	SB.SetScalarParams(ScalarNames, ScalarValues);
	SB.AllocHits();
	SB.Search();
	SB.Bench();
	//SB.WriteHits(opt(output));
	}

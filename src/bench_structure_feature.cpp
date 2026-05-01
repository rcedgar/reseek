#include "myutils.h"
#include "flat_chain.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "fastbench.h"
#include "cigar.h"
#include "chaq.h"
#include <deque>

static string s_feature;
static FastBench *s_FB;

void trunc_label(string &Label);

void read_flat_chains_idx_trunclabel(
	const string &fn,
	vector<flat_chain_t *> &chains,
	unordered_map<string, uint> &label2idx);

float flat_getlddt_muscle_some_floats(
	const uint32_t *posQs,
	const uint32_t LQ,
	const uint32_t *posTs,
	const uint32_t LT,
	const uint ncol,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint32_t *nr_considered_vec,
	uint32_t *nr_preserved_vec);

void path2posvecs(
	const string &path,
	uint loQ, uint LQ,
	uint loT, uint LT,
	vector<uint> &posQs,
	vector<uint> &posTs)
	{
	posQs.clear();
	posTs.clear();
	const uint colcount = uint(path.size());
	posQs.reserve(colcount);
	posTs.reserve(colcount);
	uint posQ = loQ;
	uint posT = loT;
	for (uint col = 0; col < colcount; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			assert(posQ < LQ);
			assert(posT < LT);
			posQs.push_back(posQ);
			posTs.push_back(posT);
			}
		if (c == 'M' || c == 'D')
			posQ++;
		if (c == 'M' || c == 'I')
			posT++;
		}
	}

uint path2posvecs2(
	const string &path,
	uint loQ, uint LQ,
	uint loT, uint LT,
	uint32_t *posQs,
	uint32_t *posTs)
	{
	const uint colcount = uint(path.size());
	uint posQ = loQ;
	uint posT = loT;
	uint m = 0;
	for (uint col = 0; col < colcount; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			assert(posQ < LQ);
			assert(posT < LT);
			posQs[m] = posQ;
			posTs[m] = posT;
			++m;
			}
		if (c == 'M' || c == 'D')
			posQ++;
		if (c == 'M' || c == 'I')
			posT++;
		}
	return m;
	}

static float get_score(
	const string &labelQ, const string &labelT,
	const sid_t *distmxQ, const sid_t *distmxT,
	const string &path, uint loQ, uint LQ, uint loT, uint LT)
	{
	if (s_feature == "lddt_old")
		{
		vector<uint> posQs, posTs;
		path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
		const uint ncol = uint(posQs.size());
		assert(posTs.size() == ncol);

		uint32_t *nr_considered_vec = myalloc(uint32_t, ncol);
		uint32_t *nr_preserved_vec = myalloc(uint32_t, ncol);
		float lddt = flat_getlddt_old_some_floats(
			posQs.data(), LQ,
			posTs.data(), LT,
			ncol,
			distmxQ,
			distmxT,
			nr_considered_vec,
			nr_preserved_vec);
		myfree(nr_considered_vec);
		myfree(nr_preserved_vec);

		float lddt_old = flat_getlddt_muscle_some_floats2(
			distmxQ, distmxT, LQ, LT, posQs, posTs);
		return lddt_old;
		}
	else if (s_feature == "lddt")
		{
		vector<uint> posQs, posTs;
		path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
		const uint ncol = uint(posQs.size());
		assert(posTs.size() == ncol);
		uint32_t *nr_considered_vec = myalloc(uint32_t, ncol);
		uint32_t *nr_preserved_vec = myalloc(uint32_t, ncol);
		float lddt = flat_getlddt_muscle_some_floats(
			posQs.data(), LQ,
			posTs.data(), LT,
			ncol,
			distmxQ,
			distmxT,
			nr_considered_vec,
			nr_preserved_vec);
		myfree(nr_considered_vec);
		myfree(nr_preserved_vec);
		return lddt;
		}
	else if (s_feature == "lddtpow")
		{
		vector<uint> posQs, posTs;
		path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
		const uint ncol = uint(posQs.size());
		assert(posTs.size() == ncol);
		uint32_t *nr_considered_vec = myalloc(uint32_t, ncol);
		uint32_t *nr_preserved_vec = myalloc(uint32_t, ncol);
		float lddt = flat_getlddt_muscle_some_floats(
			posQs.data(), LQ,
			posTs.data(), LT,
			ncol,
			distmxQ,
			distmxT,
			nr_considered_vec,
			nr_preserved_vec);
		myfree(nr_considered_vec);
		myfree(nr_preserved_vec);
		float maxL = max(LT, LQ) - 20.0f;
		if (maxL < 80)
			maxL = 80;
		uint nmatch = 0;
		for (auto c : path)
			if (c == 'M') ++nmatch;
		float score = lddt*nmatch*2.0f/powf(maxL, 0.5);
		return score;
		}
	else if (s_feature == "lddtx")
		{
		float L = (LQ + LT)/2.0f + 50;

		vector<uint> posQs, posTs;
		path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
		const uint ncol = uint(posQs.size());
		assert(posTs.size() == ncol);
		uint32_t *nr_considered_vec = myalloc(uint32_t, ncol);
		uint32_t *nr_preserved_vec = myalloc(uint32_t, ncol);
		float Lfactor = float(ncol)/L;

		asserta(distmxQ != 0 && distmxT != 0);
		float lddt = flat_getlddt_muscle_some_floats(
			posQs.data(), LQ,
			posTs.data(), LT,
			ncol, distmxQ, distmxT,
			nr_considered_vec,
			nr_preserved_vec);

		myfree(nr_considered_vec);
		myfree(nr_preserved_vec);

		float lddtx = flat_params::m_lddtx_w*lddt*500*Lfactor;
		return lddtx;
		}
	else if (s_feature == "dali")
		{
		float dali = flat_get_dali(path, 
				loQ, LQ,
				loT, LT,
				distmxQ, distmxT);
		return dali;
		}
	else if (s_feature == "dalix")
		{
		uint ncol = uint(path.size());
		float *colscores = myalloc(float, ncol);
		float dalix = flat_get_dalix(path, 
				loQ, LQ,
				loT, LT,
				distmxQ, distmxT,
				colscores);
		myfree(colscores);
		return dalix;
		}
	else
		Die("feature='%s'", s_feature.c_str());
	return 0;
	}

static void do_pair(
	flat_chain_t *chainQ,
	flat_chain_t *chainT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint loQ, uint LQ, uint loT, uint LT,
	const string &CIGAR)
	{
	const string labelQ = chainQ->m_label;
	const string labelT = chainT->m_label;

	uint domidxQ = s_FB->m_look->get_domidx(labelQ);
	uint domidxT = s_FB->m_look->get_domidx(labelT);

	const uint LQ2 = chainQ->get_length();
	const uint LT2 = chainT->get_length();
	asserta(LQ2 == LQ);
	asserta(LT2 == LT);

	string path;
	CIGARToPath(CIGAR, path, false);

	float score = get_score(
		labelQ, labelT,
		distmxQ, distmxT,
		path, loQ, LQ, loT, LT);

	uint k = s_FB->m_look->
		get_pair_idx_upper_triangle_with_diagonal(domidxQ, domidxT);
	s_FB->m_Scores[k] = score;
	}

static void read_hits(
	const string &fn,
	const unordered_map<string, uint> &label2chainidx,
	const vector<flat_chain_t *> &chains,
	const vector<sid_t *> distmxs,
	const vector<uint> &Ls)
	{
	string line;
	vector<string> flds;
	struct hit_task_t
		{
		uint chainidxQ;
		uint chainidxT;
		uint loQ;
		uint LQ;
		uint loT;
		uint LT;
		string CIGAR;
		};
	FILE *fin = OpenStdioFile(fn);
	uint64_t file_size = GetStdioFileSize64(fin);
	uint64_t last_file_pos = 0;
	uint nhit = 0;
	Progress("Reading hits 0%%\r");
	struct work_item_t
		{
		hit_task_t task;
		};
	deque<work_item_t> workq;
	const size_t MaxPendingTasks = 4096;
	mutex work_mutex;
	condition_variable work_cv;
	condition_variable space_cv;
	bool done_reading = false;
	const uint ThreadCount = GetRequestedThreadCount();
	const uint WorkerCount = (ThreadCount > 1 ? ThreadCount - 1 : 1);
	vector<thread> workers;
	workers.reserve(WorkerCount);
	for (uint w = 0; w < WorkerCount; ++w)
		{
		workers.push_back(thread([&]()
			{
			for (;;)
				{
				work_item_t work;
				{
				unique_lock<mutex> lock(work_mutex);
				work_cv.wait(lock, [&]()
					{
					return done_reading || !workq.empty();
					});
				if (workq.empty())
					{
					asserta(done_reading);
					return;
					}
				work = std::move(workq.front());
				workq.pop_front();
				space_cv.notify_one();
				}

				const uint chainidxQ = work.task.chainidxQ;
				const uint chainidxT = work.task.chainidxT;
				flat_chain_t *chainQ = chains[chainidxQ];
				flat_chain_t *chainT = chains[chainidxT];
				const uint LQ2 = chainQ->get_length();
				const uint LT2 = chainT->get_length();
				asserta(LQ2 == Ls[chainidxQ]);
				asserta(LQ2 == work.task.LQ);
				asserta(LT2 == Ls[chainidxT]);
				asserta(LT2 == work.task.LT);
				const sid_t *distmxQ = distmxs[chainidxQ];
				const sid_t *distmxT = distmxs[chainidxT];
				do_pair(
					chainQ, chainT, distmxQ, distmxT,
					work.task.loQ - 1, work.task.LQ,
					work.task.loT - 1, work.task.LT,
					work.task.CIGAR);
				}
			}));
		}
	for (;;)
		{
		bool ok = ReadLineStdioFile(fin, line);
		if (!ok)
			break;
		++nhit;
		uint64_t file_pos = GetStdioFilePos64(fin);
		if (file_pos - last_file_pos > 10e6)
			{
			double pct = file_pos*100.0/(file_size+1);
			Progress("Reading hits %.2f%%\r", pct);
			last_file_pos = file_pos;
			}

/***
                0                  1   2   3   4    5    6   7                                                      8
				q                  q qlo qhi  ql  tlo  thi  tl                                                  cigar
d1tdja3/d.58.18.2   d3nfka_/b.36.1.1  11  56  75    4   89  75                    6M1D5M3D2M4D1M2D3M8I6M1I4M2I2M39I7M
d1tdja3/d.58.18.2  d1tdja3/d.58.18.2   1  75  75    1   75  75                                                    75M
d1tdja3/d.58.18.2   d1r6ta1/a.16.1.3  69  75  75    4   14  75                                             4M2I2M2I1M
d1tdja3/d.58.18.2  d2nu8b2/d.142.1.4  10  70  75  106  230  75          10M34I2M19I5M1I1M6I4M3D1M3I6M4I7M3I7M7D3M4I5M
d1tdja3/d.58.18.2   d2gtlm1/b.61.7.1   5  59  75    4   71  75              2M1I1M3D10M6I1M14D1M5I5M16I8M1D2M4I4M1D2M
 d3nfka_/b.36.1.1   d3nfka_/b.36.1.1   1  92  92    1   92  92                                                    92M
d1tdja3/d.58.18.2   d2gc6a1/c.59.1.2   1  75  75    2  114  75     3M4I4M4D11M3I11M1D7M18I1M8I2M15I3M1I1M1D5M9D3M4I9M
 d3nfka_/b.36.1.1   d1r6ta1/a.16.1.3  66  77  92   26   43  92                                                 5M6I7M
d1tdja3/d.58.18.2  d2ziba_/d.169.1.1  12  63  75   19  114  75  3M1D3M1I6M1I8M3D3M17I3M11I4M1I2M2I2M1D3M11I5M4I1M1I4M
d1tdja3/d.58.18.2   d1mzka_/b.26.1.2   1  54  75    2   98  75                      2M7D4M1D8M8D2M10I7M21I5M27I4M1I6M
***/

		Split(line, flds, '\t');
		asserta(flds.size() == 9);
		string labelQ = flds[0];
		string labelT = flds[1];

		uint loQ = StrToUint(flds[2]);
		uint LQ = StrToUint(flds[4]);

		uint loT = StrToUint(flds[5]);
		uint LT = StrToUint(flds[7]);
		const string &CIGAR = flds[8];

		trunc_label(labelQ);
		trunc_label(labelT);

		unordered_map<string, uint>::const_iterator iterQ =
			label2chainidx.find(labelQ);
		unordered_map<string, uint>::const_iterator iterT =
			label2chainidx.find(labelT);

		if (iterQ == label2chainidx.end()) Die("Chain not found >%s", labelQ.c_str());
		if (iterT == label2chainidx.end()) Die("Chain not found >%s", labelT.c_str());

		asserta(loQ > 0);
		asserta(loT > 0);
		hit_task_t task;
		task.chainidxQ = iterQ->second;
		task.chainidxT = iterT->second;
		task.loQ = loQ;
		task.LQ = LQ;
		task.loT = loT;
		task.LT = LT;
		task.CIGAR = CIGAR;
		{
		unique_lock<mutex> lock(work_mutex);
		space_cv.wait(lock, [&]()
			{
			return workq.size() < MaxPendingTasks;
			});
		work_item_t work;
		work.task = task;
		workq.push_back(std::move(work));
		}
		work_cv.notify_one();
		}
	{
	unique_lock<mutex> lock(work_mutex);
	done_reading = true;
	}
	work_cv.notify_all();
	for (uint w = 0; w < WorkerCount; ++w)
		workers[w].join();
	Progress("Reading hits 100%%   \n");
	}

void cmd_bench_structure_feature()
	{
	asserta(optset_input);
	asserta(optset_feature);
	asserta(optset_lookup);
	s_feature = opt(feature);

	FastBench FB;
	FB.ReadLookup(opt(lookup));
	FB.Alloc();

	s_FB = &FB;

	vector<flat_chain_t *> chains;
	unordered_map<string, uint> label2chainidx;
	read_flat_chains_idx_trunclabel(opt(input), chains, label2chainidx);
	const uint nchain = uint(chains.size());

	vector<sid_t *> distmxs;
	vector<uint32_t> Ls;
	distmxs.reserve(nchain);
	Ls.reserve(nchain);
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		ProgressStep(chainidx, nchain, "Distance mxs");
		const flat_chain_t *chain = chains[chainidx];
		const uint L = chain->get_length();
		Ls.push_back(L);
		sid_t *distmx = myalloc(sid_t, L*flat_params::m_distmx_bandwidth);
		chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);
		distmxs.push_back(distmx);
		}

	read_hits(g_Arg1, label2chainidx, chains, distmxs, Ls);
	FB.SetScoreOrder();
	FB.Bench();
	}
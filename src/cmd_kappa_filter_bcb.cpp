#include "myutils.h"
#include "kappa_filter.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "bitdope.h"
#include "lookup.h"

static void bench(
	uint filter_secs,
	const BCAData &QDB,
	const BCAData &TDB)
	{
	if (!optset_dope) return;

	const uint QSeqCount = QDB.GetChainCount();
	const uint TSeqCount = TDB.GetChainCount();

	bitdope dope;
	lookup look;

	const string &lookupfn =
		(optset_lookup ? opt(lookup) : "../data/scop40x.lookup");
	look.from_tsv(lookupfn);
	dope.m_look = &look;
	dope.from_file(opt(dope));
	dope.set_square();
	ProgressLog("dope %s hits\n", FloatToStr(dope.m_nhit));

	uint npass = 0;
	uint nindope = 0;
	const vector<vector<uint16_t> > &QueryIdxToTopScoreVec =
		kappa_filter::m_RSB.m_QueryIdxToTopScoreVec;
	asserta(QueryIdxToTopScoreVec.size() == QSeqCount);
	for (uint qidx = 0; qidx < QSeqCount; ++qidx)
		{
		const string &q = QDB.GetLabel(qidx);
		uint qdomidx = look.get_domidx(q);
		const vector<uint16_t> &row = QueryIdxToTopScoreVec[qidx];
		for (uint tidx = 0; tidx < TSeqCount; ++tidx)
			{
			uint16_t score = row[tidx];
			if (score > 0)
				{
				const string &t = TDB.GetLabel(tidx);
				uint tdomidx = look.get_domidx(t);
				++npass;
				if (dope.in_square_ij(qidx, tidx))
					++nindope;
				}
			}
		}
	double pct = GetPct(nindope, 2*dope.m_nhit);

	Progress("pct=%.1f", pct);
	Progress(" secs=%u", filter_secs);
	Progress(" pattern=%s", flat_params::m_kappa_pattern.c_str());
	Progress(" kmer=%d", flat_params::m_kappa_min_kmerpairscore);
	Progress(" diag=%d", flat_params::m_kappa_min_diagscore);
	Progress(" npass=%u", npass);
	Progress("\n");

	Log("@FEV@");
	Log("\tpct=%.1f", pct);
	Log("\tsecs=%u", filter_secs);
	Log("\tpattern=%s", flat_params::m_kappa_pattern.c_str());
	Log("\tkmer=%d", flat_params::m_kappa_min_kmerpairscore);
	Log("\tdiag=%d", flat_params::m_kappa_min_diagscore);
	Log("\tnpass=%u", npass);
	Log("\n");
	}

void cmd_kappa_filter_bcb()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt(db);

	BCAData QBCA;
	BCAData TBCA;

	QBCA.Open(QFN);
	asserta(QBCA.m_HasNuSequences);

	uint8_t **query_kappa_codeseqs = 0;
	uint *query_lengths = 0;
	QBCA.make_kappa_codeseqs(&query_kappa_codeseqs, &query_lengths);

	TBCA.Open(DBFN);

	const uint QSeqCount = QBCA.GetChainCount();
	const uint TSeqCount = QBCA.GetChainCount();
	decide_query_or_db_kmer_neighborhood(QSeqCount, TSeqCount);

	kappa_filter::init_kappa();
	kappa_filter::m_RSB.Init(QSeqCount);

	kappa_dex QKmerIndex;
	QKmerIndex.Init();

	const uint k = flat_params::m_kappa_kmer_nrones;
	const kappa_mermx &GetKappaMerMx(uint k);
	const kappa_mermx &ScoreMx = GetKappaMerMx(k);
	asserta(ScoreMx.m_k == k);

	QKmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	QKmerIndex.m_MinKmerSelfScore =  flat_params::m_kappa_min_diagscore;
	QKmerIndex.from_codeseqs(query_kappa_codeseqs, query_lengths,
		QBCA.m_Labels, QSeqCount);
#if DEBUG
	QKmerIndex.Validate();
#endif
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == flat_params::m_kappa_dict_size);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	kappa_filter::m_ptrScoreMx = &ScoreMx;
	kappa_filter::m_ptrQKmerIndex = &QKmerIndex;

	kappa_seqsource db_ss;
	db_ss.OpenBCB(TBCA);

	time_t t_start = time(0);
	kappa_filter::run_filter(
		query_kappa_codeseqs, query_lengths, QSeqCount, db_ss);
	time_t t_end = time(0);
	uint filter_secs = uint(t_end - t_start);
	bench(filter_secs, QBCA, TBCA);
	}

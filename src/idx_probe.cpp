#include "myutils.h"
#include "kappa_dex.h"
#include "kappa_seqsource.h"
#include "kappa_mermx.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "sort.h"
#include "bcadata.h"

// Count-only DB-index probe: query k-mer neighborhoods vs exact kappa_dex.
//   reseek -idx_probe query.fa -db db.kdx -output probe.tsv
// TSV: query_label  target_idx  hit_count  (non-zero targets, count desc)
void cmd_idx_probe()
	{
	asserta(optset_db);
	asserta(optset_output);

	kappa_dex Index;
	Index.FromFile(opt(db));
	asserta(Index.m_nseq > 0);
	asserta(Index.m_DictSize > 0);

	const uint k = Index.m_k;
	const kappa_mermx &GetKappaMerMx(uint k);
	const kappa_mermx &ScoreMx = GetKappaMerMx(k);
	asserta(ScoreMx.m_k == k);
	asserta(ScoreMx.m_AS_pow[k] == Index.m_DictSize);

	Index.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	const int MinScore = Index.m_MinKmerSelfScore;

	uint *NeighborKmers = myalloc(uint, Index.m_DictSize);
	uint32_t *Counts = myalloc(uint32_t, Index.m_nseq);
	uint32_t *TargetIdxs = myalloc(uint32_t, Index.m_nseq);
	uint32_t *Scores_nz = myalloc(uint32_t, Index.m_nseq);
	unsigned *OrderBuf = myalloc(unsigned, Index.m_nseq);

	BCAData qbcb;
	kappa_seqsource QSS;
	uint8_t **q_codeseqs = 0;
	uint *q_lengths = 0;
	uint q_nseq = 0;
	const vector<string> *q_labels = 0;
	bool q_from_bcb = EndsWith(g_Arg1, ".bcb");
	if (q_from_bcb)
		{
		qbcb.Open(g_Arg1);
		asserta(qbcb.m_HasNuSequences);
		qbcb.make_kappa_codeseqs(&q_codeseqs, &q_lengths);
		q_nseq = qbcb.GetChainCount();
		q_labels = &qbcb.m_Labels;
		}
	else
		QSS.OpenFasta(g_Arg1);

	const char *dbname = opt(db);
	FILE *fout = CreateStdioFile(opt(output));
	fprintf(fout, "# idx_probe count-only\n");
	fprintf(fout, "# index\t%s\n", dbname);
	fprintf(fout, "# query\t%s\n", g_Arg1.c_str());
	fprintf(fout, "# k\t%u\tK\t%u\tdict\t%u\tnseq\t%u\tminself\t%d\n",
		Index.m_k, Index.m_K, Index.m_DictSize, Index.m_nseq, MinScore);
	fprintf(fout, "# fields\tquery_label\ttarget_idx\thit_count\n");

	ObjMgr OM;
	SeqInfo *SI = OM.GetSeqInfo();
	uint qidx = 0;
	uint64 TotalHits = 0;
	uint64 TotalNonZero = 0;

	auto ProbeOne = [&](const char *Label, const byte *Seq, uint QL)
		{
		if (QL < Index.m_K)
			return;

		zero_array(Counts, Index.m_nseq);
		Index.SetSeq(qidx, Label, Seq, QL);
		const vector<uint> &QKmers = Index.m_Kmers;
		const uint NK = SIZE(QKmers);

		for (uint qi = 0; qi < NK; ++qi)
			{
			uint QKmer = QKmers[qi];
			if (QKmer == UINT_MAX)
				continue;
			asserta(QKmer < Index.m_DictSize);
			if (Index.m_KmerSelfScores[QKmer] < MinScore)
				continue;

			const uint n = ScoreMx.GetHighScoringKmers(QKmer, short(MinScore),
				NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				const uint Nbr = NeighborKmers[j];
				asserta(Nbr < Index.m_DictSize);
				const uint64_t RowSize = Index.GetRowSize(Nbr);
				if (RowSize == 0)
					continue;
				uint64_t DataOffset = Index.GetRowStart(Nbr);
				for (uint64_t c = 0; c < RowSize; ++c)
					{
					uint32_t TSeqIdx;
					uint16_t TPos;
					Index.Get(DataOffset++, TSeqIdx, TPos);
					asserta(TSeqIdx < Index.m_nseq);
					Counts[TSeqIdx] += 1;
					}
				}
			}

		uint n_nz = 0;
		for (uint t = 0; t < Index.m_nseq; ++t)
			{
			if (Counts[t] == 0)
				continue;
			TargetIdxs[n_nz] = t;
			Scores_nz[n_nz] = Counts[t];
			TotalHits += Counts[t];
			++n_nz;
			}
		TotalNonZero += n_nz;

		if (n_nz > 0)
			QuickSortOrderDesc(Scores_nz, n_nz, OrderBuf);

		for (uint i = 0; i < n_nz; ++i)
			{
			const uint oi = OrderBuf[i];
			fprintf(fout, "%s\t%u\t%u\n", Label, TargetIdxs[oi], Scores_nz[oi]);
			}
		ProgressLog("query %u %s  nz=%u\n", qidx, Label, n_nz);
		};

	if (q_from_bcb)
		{
		for (uint i = 0; i < q_nseq; ++i)
			{
			ProbeOne((*q_labels)[i].c_str(), q_codeseqs[i], q_lengths[i]);
			++qidx;
			}
		qbcb.Close();
		}
	else
		{
		for (;;)
			{
			bool ok = QSS.GetNext(SI);
			if (!ok)
				break;
			ProbeOne(SI->m_Label, SI->m_Seq, SI->m_L);
			++qidx;
			}
		QSS.Close();
		}

	CloseStdioFile(fout);
	myfree(NeighborKmers);
	myfree(Counts);
	myfree(TargetIdxs);
	myfree(Scores_nz);
	myfree(OrderBuf);

	ProgressLog("idx_probe done  queries=%u  nonzero_pairs=%llu  hit_incs=%llu\n",
		qidx,
		(unsigned long long) TotalNonZero,
		(unsigned long long) TotalHits);
	}

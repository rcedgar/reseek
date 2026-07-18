#include "myutils.h"
#include "kappa_dex.h"
#include "kappa_seqsource.h"
#include "kappa_mermx.h"
#include "kappa_filter.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "twohitdiag.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "bcadata.h"
#include <set>

void set_default_stats();

// DB-index pre-HSP seed diagonals (dual of query-index -dump_prefilter_prehsp).
//   reseek -idx_prefilter_prehsp query.bcb -db db.kdx -output diags.tsv
// Optional: -input2 db.bcb|fasta for target labels; -twohitdiag / -onehitdiag
// TSV: q_label  t_label  qidx  tidx  diag
void cmd_idx_prefilter_prehsp()
	{
	asserta(optset_db);
	asserta(optset_output);

	set_default_stats();
	flat_params params;
	params.init_from_cmdline();
	flat_params::init_kappa();
	params.logme();

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
	const int MinScore = flat_params::m_kappa_min_kmerpairscore;
	if (MinScore != Index.m_MinKmerSelfScore)
		ProgressLog("Warning: -kappa_minkmerscore %d != index MinKmerSelfScore %d\n",
			MinScore, Index.m_MinKmerSelfScore);

	if (flat_params::m_kappa_onehitdiag && Index.m_nseq >= UINT16_MAX)
		Die("-onehitdiag idx_prefilter_prehsp requires db nseq < 65535");

	uint *NeighborKmers = myalloc(uint, Index.m_DictSize);

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

	vector<string> t_labels_storage;
	const vector<string> *t_labels = 0;
	BCAData tbcb;
	if (optset_input2)
		{
		const string &tfn = opt(input2);
		if (EndsWith(tfn, ".bcb"))
			{
			tbcb.Open(tfn);
			asserta(tbcb.GetChainCount() == Index.m_nseq);
			t_labels = &tbcb.m_Labels;
			}
		else
			{
			kappa_seqsource TSS;
			TSS.OpenFasta(tfn);
			ObjMgr OM;
			SeqInfo *SI = OM.GetSeqInfo();
			t_labels_storage.resize(Index.m_nseq);
			for (uint i = 0; i < Index.m_nseq; ++i)
				{
				bool ok = TSS.GetNext(SI);
				asserta(ok);
				t_labels_storage[i] = SI->m_Label;
				}
			TSS.Close();
			t_labels = &t_labels_storage;
			}
		}

	FILE *fout = CreateStdioFile(opt(output));
	kappa_filter::write_prehsp_tsv_header(fout, "db_index");
	fprintf(fout, "# index\t%s\n", opt(db));
	fprintf(fout, "# query\t%s\n", g_Arg1.c_str());
	if (optset_input2)
		fprintf(fout, "# target_labels\t%s\n", opt(input2));

	TwoHitDiag Bag;
	set<uint32_t> OneHit;
	ObjMgr OM;
	SeqInfo *SI = OM.GetSeqInfo();
	uint qidx = 0;
	uint64 TotalHits = 0;

	auto ProcessOne = [&](const char *Label, const byte *Seq, uint QL)
		{
		if (QL < Index.m_K)
			return;

		Bag.Reset();
		OneHit.clear();
		Index.SetSeq(qidx, Label, Seq, QL);
		const vector<uint> &QKmers = Index.m_Kmers;
		const uint NK = SIZE(QKmers);

		for (uint QPos = 0; QPos < NK; ++QPos)
			{
			uint QKmer = QKmers[QPos];
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
					const uint16_t Diag = uint16_t(QL + TPos - QPos - 1);
					if (Diag > m_Mask14)
						continue;
					if (flat_params::m_kappa_onehitdiag)
						{
						asserta(TSeqIdx < UINT16_MAX);
						OneHit.insert((uint32_t(TSeqIdx) << 16) | uint32_t(Diag));
						}
					else
						Bag.Add(TSeqIdx, Diag);
					}
				}
			}

		auto Emit = [&](uint tidx, uint16_t diag)
			{
			const char *tlabel = "-";
			if (t_labels != 0 && tidx < SIZE(*t_labels))
				tlabel = (*t_labels)[tidx].c_str();
			kappa_filter::write_prehsp_hit(fout, Label, tlabel, qidx, tidx, diag);
			++TotalHits;
			};

		if (flat_params::m_kappa_onehitdiag)
			{
			for (set<uint32_t>::const_iterator iter = OneHit.begin();
				 iter != OneHit.end(); ++iter)
				{
				const uint32_t pair = *iter;
				Emit(pair >> 16, uint16_t(pair & 0xffff));
				}
			}
		else
			{
			if (flat_params::m_kappa_twohitdiag)
				Bag.SetDupes();
			else
				Bag.SetUniqueFine();
			for (uint i = 0; i < Bag.m_DupeCount; ++i)
				Emit(Bag.m_DupeSeqIdxs[i], Bag.m_DupeDiags[i]);
			}
		};

	if (q_from_bcb)
		{
		for (uint i = 0; i < q_nseq; ++i)
			{
			ProcessOne((*q_labels)[i].c_str(), q_codeseqs[i], q_lengths[i]);
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
			ProcessOne(SI->m_Label, SI->m_Seq, SI->m_L);
			++qidx;
			}
		QSS.Close();
		}

	if (optset_input2 && EndsWith(string(opt(input2)), ".bcb"))
		tbcb.Close();

	CloseStdioFile(fout);
	myfree(NeighborKmers);

	ProgressLog("idx_prefilter_prehsp done  queries=%u  prehsp_hits=%llu\n",
		qidx, (unsigned long long) TotalHits);
	ProgressLog("Compare to: flat_search_kappa -idxt -dump_prefilter_prehsp ... (same -twohitdiag/-onehitdiag)\n");
	}

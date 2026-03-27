#include "myutils.h"
#include "seqdb.h"
#include "flat_chain.h"
#include "lookup.h"
#include "alpha.h"
#include "entropy.h"

void entropy::parse_fa2(
	SeqDB &DB,
	vector<uint> profidxqs,
	vector<uint> profidxts,
	vector<vector<uint> > &posvecq,
	vector<vector<uint> > &posvect)
	{
	profidxqs.clear();
	profidxts.clear();
	posvecq.clear();
	posvect.clear();

	DB.SetLabelToIndex();
	const uint nseq = DB.GetSeqCount();
	asserta(nseq%2 == 0);
	const uint npair = nseq/2;

	profidxqs.reserve(npair);
	profidxts.reserve(npair);
	posvecq.reserve(npair);
	posvect.reserve(npair);

	vector<string> flds;
	for (uint pairidx = 0; pairidx < npair; ++pairidx)
		{
		uint seqidxq = 2*pairidx;
		uint seqidxt = seqidxq + 1;

		string labelq = DB.GetLabel(seqidxq);
		const string &labelt = DB.GetLabel(seqidxt);

		Split(labelq, flds, ' ');
		labelq = flds[0];

		map<string, uint>::const_iterator iterq = m_label2idx.find(labelq);
		map<string, uint>::const_iterator itert = m_label2idx.find(labelt);
		asserta(iterq != m_label2idx.end());
		asserta(itert != m_label2idx.end());
		uint profidxq = iterq->second;
		uint profidxt = itert->second;

		const string &rowq = DB.GetSeq(seqidxq);
		const string &rowt = DB.GetSeq(seqidxt);
		const uint ncol = uint(rowq.size());
		asserta(rowt.size() == ncol);

		vector<uint> posqs;
		vector<uint> posts;
		posqs.reserve(ncol);
		posts.reserve(ncol);
		uint posq = 0;
		uint post = 0;
		for (uint col = 0; col < ncol; ++col)
			{
			char q = rowq[col];
			char t = rowq[col];
			if (isupper(q) && isupper(t))
				{
				posqs.push_back(posq);
				posts.push_back(post);
				}
			if (!isgap(q))
				++posq;
			if (!isgap(t))
				++post;
			}

		profidxqs.push_back(profidxq);
		profidxts.push_back(profidxt);
		posvecq.push_back(posqs);
		posvect.push_back(posts);
		}
	}

void entropy::load_fa2s(
	const string &lookupfn,
	const string &tpfa2fn,
	const string &fpfa2fn)
	{
	m_TPDB.FromFasta(opt(fasta2_tp), true);
	m_FPDB.FromFasta(opt(fasta2_fp), true);

	parse_fa2(m_TPDB, 
		m_tp_profidxqs, m_tp_profidxts,
		m_tp_posvecq, m_tp_posvect);

	parse_fa2(m_FPDB,
		m_fp_profidxqs, m_fp_profidxts,
		m_fp_posvecq, m_fp_posvect);
	}

void cmd_train_fa2auc()
	{
	asserta(optset_lookup);
	asserta(optset_fasta2_tp);
	asserta(optset_fasta2_fp);

	entropy E;

	const string &filesfn = g_Arg1;
	vector<string> lines;
	ReadLinesFromFile(filesfn, lines);

	const size_t nfeat = lines.size();
	vector<string> fafns;
	vector<string> logoddsfns;
	vector<string> flds;
	for (size_t fi = 0; fi < nfeat; ++fi)
		{
		Split(lines[fi], flds, '\t');
		asserta(flds.size() == 2);
		fafns.push_back(flds[0]);
		logoddsfns.push_back(flds[1]);
		}

	E.load_profiles(fafns);
	E.read_logoddsvec(logoddsfns);
	E.load_fa2s(opt(lookup), opt(fasta2_tp), opt(fasta2_fp));

	Progress("done.\n");
	}

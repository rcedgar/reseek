#include "myutils.h"
#include "seqdb.h"
#include "flat_helpers.h"

// Make FASTA for letters in a compond e.g. aa4+pm2
// with combined alphabet size <= 36
void cmd_flat_merge_features()
	{
	const string &feature_names_str = g_Arg1;
	vector<string> feature_names;
	Split(feature_names_str, feature_names, '+');
	const uint nfeat = uint(feature_names.size());

	asserta(optset_fasta);
	const string &fastafn = opt(fasta);
	FILE *ffa = CreateStdioFile(fastafn);
	const string &fafnpattern = opt(fapattern);
	const string &logoddsfnpattern = opt(mxpattern);

	vector<string> fafns(nfeat);
	vector<uint> alpha_sizes(nfeat);
	vector<uint> factors;
	uint compound_alpha_size = 1;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const string &feature_name = feature_names[fi];
		uint alpha_size = get_alpha_size_from_feature_name(feature_name);
		asserta(alpha_size != 20); // aa is special case for chartoletter
		factors.push_back(compound_alpha_size);
		compound_alpha_size *= alpha_size;
		asserta(compound_alpha_size <= 36); // uses g_CharToLetterMu
		alpha_sizes[fi] = alpha_size;

		make_fn_pattern(
			fafnpattern,
			feature_name,
			fafns[fi]);
		}

	vector<SeqDB *> DBs;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		SeqDB *DB = new SeqDB;
		DB->FromFasta(fafns[fi]);
		DB->SetLabelToIndex();
		DBs.push_back(DB);
		}
	const map<string, uint> &label2idx = DBs[0]->m_LabelToIndex;
	vector<string> seqs(nfeat);
	for (map<string, uint>::const_iterator iter = label2idx.begin();
		iter != label2idx.end(); ++iter)
		{
		const string &label = iter->first;
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			DBs[fi]->GetSeqByLabel(label, seqs[fi]);
			if (fi > 0)
				asserta(seqs[fi].size() == seqs[0].size());
			}
		const size_t L = seqs[0].size();
		string compound_seq;
		for (uint pos = 0; pos < L; ++pos)
			{
			uint compound_letter = 0;
			for (uint fi = 0; fi < nfeat; ++fi)
				{
				char c = seqs[fi][pos];
				uint letter = g_CharToLetterMu[c];
				asserta(letter < alpha_sizes[fi]);
				compound_letter += letter*factors[fi];
				}
			asserta(compound_letter < compound_alpha_size);
			compound_seq.push_back(g_LetterToCharMu[compound_letter]);
			}
		SeqToFasta(ffa, label, compound_seq);
		}
	CloseStdioFile(ffa);
	}

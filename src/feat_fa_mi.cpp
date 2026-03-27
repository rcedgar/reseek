#include "myutils.h"
#include "seqdb.h"
#include "chaq.h"

void cmd_feat_fa_mi()
	{
	asserta(optset_alpha_size1);
	asserta(optset_alpha_size2);
	const uint alpha_size1 = opt(alpha_size1);
	const uint alpha_size2 = opt(alpha_size2);
	const uint8_t *char2letter1 = chaq::get_char2letter(alpha_size1);
	const uint8_t *char2letter2 = chaq::get_char2letter(alpha_size2);

	SeqDB DB1, DB2;
	DB1.FromFasta(g_Arg1);
	DB2.FromFasta(opt(input2));

	DB2.SetLabelToIndex();

	vector<uint> counts1(alpha_size1);
	vector<uint> counts2(alpha_size2);
	vector<vector<uint> > joint_counts(alpha_size1);
	for (uint code1 = 0; code1 < alpha_size1; ++code1)
		joint_counts[code1].resize(alpha_size2);

	const uint nseq = DB1.GetSeqCount();
	uint N = 0;
	uint missing = 0;
	for (uint idx1 = 0; idx1 < nseq; ++idx1)
		{
		const string &label1 = DB1.GetLabel(idx1);
		uint idx2 = DB2.GetSeqIndex(label1);
		if (idx2 == UINT_MAX)
			{
			++missing;
			continue;
			}
		const string &seq1 = DB1.GetSeq(idx1);
		const string &seq2 = DB2.GetSeq(idx2);
		const uint L = DB1.GetSeqLength(idx1);
		asserta(DB2.GetSeqLength(idx2) == L);

		for (uint pos = 0; pos < L; ++pos)
			{
			char c1 = seq1[pos];
			char c2 = seq2[pos];
			uint8_t code1 = char2letter1[c1];
			uint8_t code2 = char2letter2[c2];
			asserta(code1 < alpha_size1);
			asserta(code2 < alpha_size2);

			counts1[code1] += 1;
			counts2[code2] += 1;
			joint_counts[code1][code2] += 1;
			++N;
			}
		}
	ProgressLog("%u missing, N=%s\n", missing, IntToStr(N));

	vector<double> freqs1(alpha_size1);
	vector<double> freqs2(alpha_size2);
	vector<vector<double> > joint_freqs(alpha_size1);
	for (uint code1 = 0; code1 < alpha_size1; ++code1)
		joint_freqs[code1].resize(alpha_size2);

	double sum = 0;
	for (uint i = 0; i < alpha_size1; ++i)
		{
		double freq = double(counts1[i])/N;
		sum += freq;
		freqs1[i] = freq;
		}
	asserta(sum > 0.99 && sum < 1.01);

	sum = 0;
	for (uint i = 0; i < alpha_size2; ++i)
		{
		double freq = double(counts2[i])/N;
		sum += freq;
		freqs2[i] = freq;
		}
	asserta(sum > 0.99 && sum < 1.01);

	sum = 0;
	for (uint i = 0; i < alpha_size1; ++i)
		{
		for (uint j = 0; j < alpha_size2; ++j)
			{
			double freq = double(joint_counts[i][j])/N;
			sum += freq;
			joint_freqs[i][j] = freq;
			}
		}
	asserta(sum > 0.99 && sum < 1.01);

	double MI = 0;
	for (uint i = 0; i < alpha_size1; ++i)
		{
		double freq1 = freqs1[i];
		for (uint j = 0; j < alpha_size2; ++j)
			{
			double joint_freq = joint_freqs[i][j];
			double freq2 = freqs2[j];

			if (joint_freq > 1e-6 && freq1 > 1e-6 && freq2 > 1e-6)
				{
				MI += joint_freq*log(joint_freq/(freq1*freq2));
				asserta(!isnan(MI));
				}
			}
		}
	ProgressLog("MI=%.3g feat1=%s feat2=%s\n", MI, g_Arg1.c_str(), opt(input2));
	}

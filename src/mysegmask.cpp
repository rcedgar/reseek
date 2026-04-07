#include "myutils.h"
#include "seqdb.h"
#include "flat_helpers.h"
#include "quarts.h"

/***
segmasker defaults:
    Window Size (W): 12 (amino acids/nucleotides)
    Low Cutoff (K1​): 2.2 bits
    High Cutoff (K2​): 2.5 bits

Triggering: It scans the sequence with the window of size W.
If the Shannon entropy of a window is less than or equal to the
low cutoff (K1​), a "low-complexity" seed is triggered.

Extension: Once triggered, the segment is extended in both
directions as long as the entropy of the resulting merged
segment remains below the high cutoff (K2​).
***/

static double get_entropy(
	const string &seq,
	uint startpos,
	uint w)
	{
	unordered_map<char, uint> c2n;
	assert(startpos + w <= seq.size());
	for (uint pos = startpos; pos < startpos+w; ++pos)
		{
		char c = seq[pos];
		unordered_map<char, uint>::iterator iter = c2n.find(c);
		if (iter == c2n.end())
			c2n[c] = 1;
		else
			c2n[c] += 1;
		}

	double H = 0;
	for (unordered_map<char, uint>::iterator iter = c2n.begin();
		iter != c2n.end(); ++iter)
		{
		uint n = iter->second;
		double p = double(n)/w;
		H += -p*log(p);
		}
	return H;
	}

void cmd_mysegmask()
	{
	const string &fastafn = g_Arg1;
	const uint alpha_size = optset_alpha_size ? opt(alpha_size) : 20;
	FILE *f = CreateStdioFile(opt(output));

	const uint W = 12;
	string seq;
	for (uint i = 0; i < W; ++i) seq += 'A' + i;
	const double maxH = get_entropy(seq, 0, W);
	ProgressLog("maxH %.4f\n", maxH);

	SeqDB DB;
	DB.FromFasta(fastafn);
	const uint8_t *char2letter = get_char2letter(alpha_size);

	uint nmasked = 0;
	uint N = 0;
	const uint nseq = DB.GetSeqCount();
	vector<float> Hs;
	for (uint seqidx = 0; seqidx < nseq; ++seqidx)
		{
		double pct = (nmasked*100.0)/(N+1);
		ProgressStep(seqidx, nseq, "Masking %u (%.3g%%)", nmasked, pct);
		const string &seq = DB.GetSeq(seqidx);
		const uint L = uint(seq.size());
		N += L;
		for (uint startpos = 0; startpos + W <= L; ++startpos)
			{
			double H = get_entropy(seq, startpos, W);
			Hs.push_back(float(H));
			if (H < 1.42)
				++nmasked;
			}
		}

	double pct = (nmasked*100.0)/N;
	ProgressLog("Masking %u seqs, %u (%.3g%%)",
		nseq, nmasked, pct);

	QuartsFloat QF;
	GetQuartsFloat(Hs, QF);
	QF.ProgressLogMe();
	}
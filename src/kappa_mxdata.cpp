#include "myutils.h"
#include "kappa_mermx.h"
#include "flat_params.h"

/***
Reduced alphabet sec4
---------------------
	C:/src/reseek_tune2/kappa_prefiter_try_sec4_specs/report.txt
	Tied for best (small differences only):
		prefilter_kappa_AB-CDGINZ-EFJLKMORSTacdef-HPQUVWXYb.tm_0.8_1.0.log:@FEV@
			pct=69.7	secs=11	pattern=11010001	kmer=50	diag=100
																 ^^^
										     Note diag = 2*kmer, ~no
									    effect might as well be zero

	Best parameters
		sec4_groups		AB-CDGINZ-EFJLKMORSTacdef-HPQUVWXYb
		train fa2		tm_0.8_1.0

Train float logodds
-------------------
	C:/src/reseek_tune2/bash/sweep_kappa_params_setup.bash
	/src/reseek_tune2 [a407b21] (HEAD -> main, upstream/main) sweep_kappa_params
	reseek [81d6520]

	fa2=tm_0.8_1.0

	reseek \
		-kappa_fasta ../data/scop40x.bca \
		-output kappa.fasta \
		-log kappa.fasta.log

	reseek \
		-flat_train_discrete kappa.fasta \
		-fasta2_tp ../big_fa2/tp.$fa2.fa2 \
		-alpha_size 32 \
		-log train_kappa.log \
		-background_style unaln \
		-output kappa.logodds

	=> kappa.logodds (float):
		logodds  32
		3.757  2.422  2.355  3.294  0.2406  -0.5578  -0.7259  ...

	-scalef 5 (convert to integer using C++ round() equiv. python round()) =>
		19  12  12  16  1  -3  -4 ...
***/

int16_t kappa32_flat_logodds[1024] = {
  19,  12,  12,  16,   1,  -3,  -4,  -2, -30, -24, -30,   0, -35, -36, -48,   0,  -2,  -5,  -3,  -7, -17, -18, -17, -19, -17, -17, -18, -20, -29, -30, -32, -28,
  12,  18,   7,  12,  -3,   2,  -8,  -4, -26, -20, -31,   0, -34, -39, -42,   0,  -4,  -1,  -6,  -8, -17, -14, -20, -21, -17, -15, -20, -22, -28, -25, -28, -28,
  12,   7,  18,   8,  -5,  -8,   2,  -9, -29, -31, -27,   0, -38, -43, -36,   0,  -7,  -9,   2, -14, -20, -23, -13, -25, -20, -18, -17, -23, -30, -33, -29, -32,
  16,  12,   8,  31,  -1,  -3,  -8,   6, -27, -28, -31,   0,   0,   0, -31,   0,  -1,  -3,  -4,   6, -12, -17, -18, -11, -17, -13, -16, -16, -25, -28, -29, -20,
   1,  -3,  -5,  -1,  18,  11,  11,  16, -45, -37, -46,   0, -26, -28, -33,   0, -18, -18, -17, -21,  -2,  -5,  -3,  -8, -30, -26, -33, -31, -19, -19, -19, -23,
  -3,   2,  -8,  -3,  11,  17,   7,  11, -36, -32, -48,   0, -27, -24, -30, -28, -16, -14, -19, -20,  -4,   0,  -5,  -9, -29, -26, -30, -28, -18, -15, -22, -24,
  -4,  -8,   2,  -8,  11,   7,  18,   7, -36, -37, -36,   0, -26, -28, -27,   0, -20, -23, -12, -25,  -7,  -9,   2, -13, -29, -31, -26, -34, -18, -19, -16, -22,
  -2,  -4,  -9,   6,  16,  11,   7,  28,   0,   0,   0,   0, -24, -23, -25,   0, -17, -16, -18,  -7,  -3,  -3,  -4,   5, -32, -22, -21,   0, -17, -14, -19,   0,
 -30, -26, -29, -27, -45, -36, -36,   0,  19,  13,  12,  15,   2,  -4,  -4,   0, -13, -15, -12, -14, -26, -28, -27, -30,   5,   1,   1,   1, -12, -14, -13, -14,
 -24, -20, -31, -28, -37, -32, -37,   0,  13,  21,   8,   9,  -4,   4,  -8,  -2, -13, -11, -14, -13, -27, -25, -25, -35,   1,   6,  -1,  -2, -15,  -9, -16, -12,
 -30, -31, -27, -31, -46, -48, -36,   0,  12,   8,  19,   6,  -3,  -7,   3,  -8, -17, -18, -11, -20, -34, -32, -27, -30,  -2,  -5,   5,  -8, -16, -18,  -9, -21,
   0,   0,   0,   0,   0,   0,   0,   0,  15,   9,   6,  32,  -1, -10,  -6,  13, -17, -15, -12,  -7, -29, -29, -35, -27,   1,   0,  -2,  16, -15, -15, -14,  -2,
 -35, -34, -38,   0, -26, -27, -26, -24,   2,  -4,  -3,  -1,  19,  13,  12,  15, -28, -25, -24, -23, -13, -15, -14, -15, -11, -16, -12, -15,   4,   1,   0,  -1,
 -36, -39, -43,   0, -28, -24, -28, -23,  -4,   4,  -7, -10,  13,  21,   9,   9, -25, -24, -26, -24, -15, -12, -16, -16, -13, -11, -15, -16,   1,   6,  -1,  -2,
 -48, -42, -36, -31, -33, -30, -27, -25,  -4,  -8,   3,  -6,  12,   9,  19,   8, -28, -32, -26, -29, -17, -18, -13, -19, -16, -19,  -9, -20,  -2,  -5,   6,  -7,
   0,   0,   0,   0,   0, -28,   0,   0,   0,  -2,  -8,  13,  15,   9,   8,  33, -22, -27, -29, -21, -15, -14, -17, -10, -12, -15, -16,  -4,   1,  -1,  -4,  15,
  -2,  -4,  -7,  -1, -18, -16, -20, -17, -13, -13, -17, -17, -28, -25, -28, -22,  13,   8,   7,   6,  -4,  -7,  -9,  -9,   2,  -1,  -4,  -2, -14, -15, -18, -16,
  -5,  -1,  -9,  -3, -18, -14, -23, -16, -15, -11, -18, -15, -25, -24, -32, -27,   8,  12,   4,   6,  -8,  -3, -12,  -9,  -2,   2,  -6,  -3, -16, -13, -19, -16,
  -3,  -6,   2,  -4, -17, -19, -12, -18, -12, -14, -11, -12, -24, -26, -26, -29,   7,   4,  17,   1,  -9, -11,  -1, -13,  -2,  -4,   5,  -5, -17, -17, -10, -22,
  -7,  -8, -14,   6, -21, -20, -25,  -7, -14, -13, -20,  -7, -23, -24, -29, -21,   6,   6,   1,  21,  -9,  -8, -13,   4,  -4,  -2,  -7,   7, -16, -17, -19,  -9,
 -17, -17, -20, -12,  -2,  -4,  -7,  -3, -26, -27, -34, -29, -13, -15, -17, -15,  -4,  -8,  -9,  -9,  13,   8,   7,   6, -14, -16, -17, -17,   2,  -1,  -3,  -1,
 -18, -14, -23, -17,  -5,   0,  -9,  -3, -28, -25, -32, -29, -15, -12, -18, -14,  -7,  -3, -11,  -8,   8,  13,   4,   7, -17, -14, -21, -16,  -2,   1,  -6,  -4,
 -17, -20, -13, -18,  -3,  -5,   2,  -4, -27, -25, -27, -35, -14, -16, -13, -17,  -9, -12,  -1, -13,   7,   4,  16,   1, -18, -18, -12, -23,  -2,  -4,   4,  -5,
 -19, -21, -25, -11,  -8,  -9, -13,   5, -30, -35, -30, -27, -15, -16, -19, -10,  -9,  -9, -13,   4,   6,   7,   1,  21, -18, -18, -20, -10,  -3,  -4,  -8,   8,
 -17, -17, -20, -17, -30, -29, -29, -32,   5,   1,  -2,   1, -11, -13, -16, -12,   2,  -2,  -2,  -4, -14, -17, -18, -18,  17,  11,   9,  10,   0,  -5,  -6,  -5,
 -17, -15, -18, -13, -26, -26, -31, -22,   1,   6,  -5,   0, -16, -11, -19, -15,  -1,   2,  -4,  -2, -16, -14, -18, -18,  11,  17,   7,   8,  -5,   1,  -9,  -6,
 -18, -20, -17, -16, -33, -30, -26, -21,   1,  -1,   5,  -2, -12, -15,  -9, -16,  -4,  -6,   5,  -7, -17, -21, -12, -20,   9,   7,  19,   5,  -6,  -8,   3,  -9,
 -20, -22, -23, -16, -31, -28, -34,   0,   1,  -2,  -8,  16, -15, -16, -20,  -4,  -2,  -3,  -5,   7, -17, -16, -23, -10,  10,   8,   5,  28,  -6,  -8,  -8,  10,
 -29, -28, -30, -25, -19, -18, -18, -17, -12, -15, -16, -15,   4,   1,  -2,   1, -14, -16, -17, -16,   2,  -2,  -2,  -3,   0,  -5,  -6,  -6,  16,  11,  10,  10,
 -30, -25, -33, -28, -19, -15, -19, -14, -14,  -9, -18, -15,   1,   6,  -5,  -1, -15, -13, -17, -17,  -1,   1,  -4,  -4,  -5,   1,  -8,  -8,  11,  17,   7,   8,
 -32, -28, -29, -29, -19, -22, -16, -19, -13, -16,  -9, -14,   0,  -1,   6,  -4, -18, -19, -10, -19,  -3,  -6,   4,  -8,  -6,  -9,   3,  -8,  10,   7,  19,   5,
 -28, -28, -32, -20, -23, -24, -22,   0, -14, -12, -21,  -2,  -1,  -2,  -7,  15, -16, -16, -22,  -9,  -1,  -4,  -5,   8,  -5,  -6,  -9,  10,  10,   8,   5,  27
};


static kappa_mermx *s_ptrkappaMerMx = 0;

const kappa_mermx &GetKappaMerMx(uint k)
	{
	if (s_ptrkappaMerMx != 0)
		return *s_ptrkappaMerMx;
	s_ptrkappaMerMx = new kappa_mermx;

	const uint AS = 32;

	short **MxPtrs = myalloc(short *, AS);
	for (uint i = 0; i < AS; ++i)
		{
		short *Row = myalloc(short, AS);
		for (uint j = 0; j < AS; ++j)
			Row[j] = kappa32_flat_logodds[i*AS + j];
		MxPtrs[i] = Row;
		}
	(*s_ptrkappaMerMx).Init(MxPtrs, k, AS, 2);
	return *s_ptrkappaMerMx;
	}

//	static uint read_logodds(const string &fn, vector<float> &logodds);
void load_kappa_integer_logodds(const string &fn, double scalef)
	{
	vector<float> logodds;
	if (fn == "")
		{
		extern const vector<string> g_alpha_collect_lines;
		collect C; // TODO should be one global collect object for defaults
		C.from_lines(g_alpha_collect_lines);

		const vector<string> &logodds_lines = C.get_lines("kappa32.logodds");
		uint alpha_size = flat_params::lines2logoddsmx(logodds_lines, logodds);
		asserta(alpha_size == 32);

		const uint n = alpha_size*alpha_size;
		for (uint k = 0; k < n; ++k)
			{
			const float score = logodds[k];
			asserta(score >= MIN_SANE_SCORE && score <= MAX_SANE_SCORE);
			logodds[k] = score;
			}
		}
	else
		{
		flat_params::read_logodds(fn, logodds);
		asserta(logodds.size() == 32*32);
		}
	Log("// before scale=%.1f\n", scalef);
	Log("int16_t logodds[1024] = {\n");
	for (uint i = 0; i < 32*32; ++i)
		{
		if (i > 0)
			Log(",");
		if (i > 0 && i%32 == 0)
			Log("\n");
		Log(" %.1f", logodds[i]);
		}
	Log("\n};\n");

	for (uint i = 0; i < 32*32; ++i)
		{
		double score = logodds[i]*scalef;
		int intscore = int(round(score));
		int16_t intscore16 = int16_t(intscore);
		asserta(intscore16 == intscore);
		kappa32_flat_logodds[i] = intscore16;
		}

	Log("int16_t kappa32_flat_logodds[1024] = {\n");
	for (uint i = 0; i < 32*32; ++i)
		{
		if (i > 0)
			Log(",");
		if (i > 0 && i%32 == 0)
			Log("\n");
		Log(" %3d", kappa32_flat_logodds[i]);
		}
	Log("\n};\n");
	}

void cmd_load_kappa_logodds()
	{
	asserta(optset_scalef);
	load_kappa_integer_logodds(g_Arg1, opt(scalef));
	}
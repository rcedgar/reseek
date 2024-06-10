#include "myutils.h"
#include "tma.h"
#include "seqdb.h"

void cmd_tmscore()
	{
	const string &QueryFileName = g_Arg1;
	const string &RefFileName = opt_ref;
	const string &AlnFileName = opt_input;

	SeqDB Aln;
	Aln.FromFasta(AlnFileName, true);
	asserta(Aln.IsAligned());
	Aln.SetLabelToIndex();

	TMA T;

	PDBChain Q;
	PDBChain R;
	Q.FromCal(QueryFileName);
	R.FromCal(RefFileName);

	const string &QLabel = Q.m_Label;
	const string &RLabel = R.m_Label;

	string RowQ;
	string RowR;
	Aln.GetSeqByLabel(QLabel, RowQ);
	Aln.GetSeqByLabel(RLabel, RowR);
	asserta(SIZE(RowQ) == SIZE(RowR));

	string QAcc;
	string RAcc;
	Q.GetAcc(QAcc);
	R.GetAcc(RAcc);

	T.AlignChains(Q, R);
	ProgressLog("Align: TM1=%.4f TM2=%.4f\n", T.m_TM1, T.m_TM2);

	T.CalcTMScore(Q, R, RowQ, RowR);
	ProgressLog("Score: TM1=%.4f TM2=%.4f\n", T.m_TM1, T.m_TM2);
	}

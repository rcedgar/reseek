#include "myutils.h"
#include "seqdb.h"
#include "chainreader2.h"
#include "alpha.h"
#include <set>

double GetDALIScore(const PDBChain& Q, const PDBChain& T,
	const vector<uint>& PosQs, const vector<uint>& PosTs);
void GetAlignedPositions(const string& RowQ, const string& RowR,
	vector<uint>& PosQs, vector<uint>& PosRs, vector<bool>* ptrCore);
void GetUngappedSeq(const string& Row, string& Seq);
double GetDALIZFromScoreAndLengths(double DALIScore, uint QL, uint TL);

/***
dali2tsv.py
	s = Query # 0
	s += "\t" + Chains[i] # 1
	s += "\t%.1f" % Z # 2
	s += "\t%.1f" % RMSDs[i] # 3
	s += "\t%d" % Lalis[i] # 4
	s += "\t%d" % NRs[i] # 5
	s += "\t%d" % int(PctIds[i]) # 6
	s += "\t" + QRow # 7
	s += "\t" + SRow # 8
	print(s)
***/

void cmd_daliscore_tsv()
	{
	const string &TsvFN = g_Arg1;

	PDBFileScanner FS;
	FS.Open(opt_input);
	ChainReader2 CR;
	CR.Open(FS);

	vector<PDBChain*> Chains;
	map<string, uint> LabelToChainIdx;
	for (;;)
		{
		PDBChain* ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		uint Idx = SIZE(Chains);
		Chains.push_back(ptrChain);
		LabelToChainIdx[ptrChain->m_Label] = Idx;
		}

	FILE *fIn = OpenStdioFile(TsvFN);
	FILE* fOut = CreateStdioFile(opt_output);

	string Line;
	vector<string> Fields;

	set<string> NotFound;
	uint PairCount = 0;
	double SumZ = 0;
	double SumZCore = 0;
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 9);
		const string &LabelQ = Fields[0];
		const string &LabelR = Fields[1];
		double Zin = StrToFloat(Fields[2]);
		uint Lali = StrToUint(Fields[4]);
		const string &RowQ = Fields[7];
		const string &RowR = Fields[8];

		map<string, uint>::const_iterator iterQ = LabelToChainIdx.find(LabelQ);
		map<string, uint>::const_iterator iterR = LabelToChainIdx.find(LabelR);
		asserta(iterQ != LabelToChainIdx.end());
		asserta(iterR != LabelToChainIdx.end());
		uint ChainIdxQ = iterQ->second;
		uint ChainIdxR = iterR->second;
		PDBChain &ChainQ = *Chains[ChainIdxQ];
		PDBChain &ChainR = *Chains[ChainIdxR];
		asserta(ChainQ.m_Label == LabelQ);
		asserta(ChainR.m_Label == LabelR);
		uint LQ = ChainQ.GetSeqLength();
		uint LR = ChainR.GetSeqLength();

		vector<uint> PosQs;
		vector<uint> PosRs;
		GetAlignedPositions(RowQ, RowR, PosQs, PosRs, 0);

		double Score = GetDALIScore(ChainQ, ChainR, PosQs, PosRs);
		double Z = GetDALIZFromScoreAndLengths(Score, LQ, LR);
		Log("%.1f %.1f %s %s\n", Zin, Z, LabelQ.c_str(), LabelR.c_str());
		}
	CloseStdioFile(fIn);
	CloseStdioFile(fOut);
	}

#include "myutils.h"
#include "seqdb.h"
#include "twohitdiag.h"
#include "mudex.h"
#include "mermx.h"
#include "diaghsp.h"
#include "diag.h"
#include "alpha.h"
#include "scop40bench.h"
#include "dssaligner.h"

void StrToMuLetters(const string &StrSeq, byte *Letters);

static vector<PDBChain *> m_Chains;
static vector<vector<byte> > m_MuLettersVec;
static uint m_NH2, m_NT2, m_NT3;
static SeqDB *m_Input;
static vector<vector<vector<byte> > *> m_Profiles;
static uint m_SeqCount;
static byte **m_MuSeqs = myalloc(byte *, m_SeqCount);
static uint m_DictSize;
static MuDex m_MD;
static MerMx m_MM;
static uint m_k;
static FILE *m_fOut;
static int m_MinDiagScore;
static int m_HSMinScore;

static void DoQ(uint SeqIdxQ, DSSAligner &DA, TwoHitDiag &TH, DiagHSP &DH,
				uint *HSKmers, int *SeqIdxTToBestDiagScore, int *SeqIdxTToBestDiag,
				uint *SeqIdxTs)
	{
	double FractDone = double(SeqIdxQ)/m_SeqCount;
	double IdealNT = FractDone*173058;
	double PctIdeal = m_NT3*100.0/(IdealNT + 1);
	ProgressStep(SeqIdxQ, m_SeqCount, "Searching %u %u (%.1f, %.1f%%)",
					m_NT2, m_NT3, IdealNT, PctIdeal);
	const PDBChain &ChainQ = *m_Chains[SeqIdxQ];
	const vector<vector<byte> > *ptrProfileQ = m_Profiles[SeqIdxQ];
	vector<byte> *ptrMuLettersQ = &m_MuLettersVec[SeqIdxQ];
	DA.SetQuery(ChainQ, ptrProfileQ, 0, ptrMuLettersQ, 0);

	const char *SeqQ = m_Input->GetSeq(SeqIdxQ).c_str();
	const string &LabelQ = m_Input->GetLabel(SeqIdxQ);
	uint LQ32 = m_Input->GetSeqLength(SeqIdxQ);
	asserta(LQ32 < UINT16_MAX);
	uint16_t LQ = uint16_t(LQ32);

// Initialize k-mer scan of Query
	uint Kmer = 0;
	for (uint SeqPosQ = 0; SeqPosQ < m_k-1; ++SeqPosQ)
		{
		byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		}

// m_k-mer scan of Query
	for (uint SeqPosQ = m_k-1; SeqPosQ < LQ; ++SeqPosQ)
		{
		byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		Kmer %= m_DictSize;
	// Construct high-scoring neighborhood of current m_k-mer (Kmer)
		const uint HSKmerCount =
			m_MM.GetHighScoring5mers(Kmer, m_HSMinScore, HSKmers);

		for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
			{
			uint HSKmer = HSKmers[HSKmerIdx];
		// Look up HSKmer in MuDex (Mu m_k-mer index of db)
			uint RowSize = m_MD.GetRowSize(HSKmer);
			uint DataOffset = m_MD.GetRowStart(HSKmer);
			for (uint ColIdx = 0; ColIdx < RowSize; ++ColIdx)
				{
				uint32_t SeqIdxT;
				uint16_t SeqPosT;
				m_MD.Get(DataOffset++, SeqIdxT, SeqPosT);
				if (SeqIdxT == SeqIdxQ)
					continue;
				uint LT32 = m_Input->GetSeqLength(SeqIdxT);
				asserta(LT32 < UINT16_MAX);
				uint16_t LT = uint16_t(LT32);
				diag dg(LQ, LT);
				uint16_t Diag = dg.getd(SeqPosQ, SeqPosT);
				TH.Add(SeqIdxT, Diag);
				}
			}
		}
	TH.ClearDupes();
	TH.SetDupes();
#if DEBUG
	TH.Validate(m_SeqCount, INT16_MAX);
#endif
	vector<pair<uint32_t, uint16_t> > SeqIdxDiagPairs;
	uint DupeCount = TH.m_DupeCount;
	const byte *MuSeqQ = m_MuSeqs[SeqIdxQ];
	DH.SetQ(MuSeqQ, LQ);
	uint HitCount = 0;
	for (uint i = 0; i < DupeCount; ++i)
		{
		uint32_t SeqIdxT = TH.m_DupeSeqIdxs[i];
		uint16_t Diag = TH.m_DupeDiags[i];
		const byte *MuSeqT = m_MuSeqs[SeqIdxT];
		const uint LT = m_Input->GetSeqLength(SeqIdxT);
		DH.SetT(MuSeqT, LT);
		int Lo, Len;
		int DiagScore = DH.Search(Diag, Lo, Len);
		if (DiagScore >= m_MinDiagScore)
			{
			int BestDiagScoreT = SeqIdxTToBestDiagScore[SeqIdxT];
			if (BestDiagScoreT == 0)
				{
				SeqIdxTs[HitCount++] = SeqIdxT;
				SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
				SeqIdxTToBestDiag[SeqIdxT] = Diag;
				//SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
				//SeqIdxTToBestDiagLen[SeqIdxT] = Len;
				}
			else if (DiagScore > SeqIdxTToBestDiagScore[SeqIdxT])
				{
				SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
				SeqIdxTToBestDiag[SeqIdxT] = Diag;
				//SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
				//SeqIdxTToBestDiagLen[SeqIdxT] = Len;
				}
			}
		}

	if (HitCount > 0)
		{
		for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
			{
			uint SeqIdxT = SeqIdxTs[HitIdx];
			const string &LabelT = m_Input->GetLabel(SeqIdxT);
			bool IsT = SCOP40Bench::IsTP_SF(LabelQ, LabelT);
			++m_NH2;
			if (IsT)
				++m_NT2;
			//int BestDiagScore = SeqIdxTToBestDiagScore[SeqIdxT];
			//int BestDiag = SeqIdxTToBestDiag[SeqIdxT];
			//int BestDiagLo = SeqIdxTToBestDiagLo[SeqIdxT];
			//int BestDiagLen = SeqIdxTToBestDiagLen[SeqIdxT];

			const PDBChain *ChainT = m_Chains[SeqIdxT];
			const vector<vector<byte> > *ptrProfileT = m_Profiles[SeqIdxT];
			vector<byte> *ptrMuLettersT = &m_MuLettersVec[SeqIdxT];
			DA.SetTarget(*ChainT, ptrProfileT, 0, ptrMuLettersT, 0);
			DA.AlignQueryTarget();
			float Evalue = DA.m_EvalueA;
			if (IsT && Evalue <= 10)
				++m_NT3;
			if (m_fOut != 0)
				fprintf(m_fOut, "%s\t%s\t%.3g\n",
						LabelQ.c_str(), LabelT.c_str(), Evalue);
			}
		}

	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		uint SeqIdxT = SeqIdxTs[HitIdx];
		SeqIdxTToBestDiagScore[SeqIdxT] = 0;
		}
#if DEBUG
	{
	for (uint SeqIdx = 0; SeqIdx < m_SeqCount; ++SeqIdx)
		{
		assert(SeqIdxTToBestDiagScore[SeqIdx] == 0);
		}
	}
#endif
	TH.Reset();
	}

void cmd_twohit3()
	{
	asserta(optset_hsminscore);
	asserta(optset_mindiagscore);
	asserta(optset_ref);

	ReadChains(opt_ref, m_Chains);
	const uint ChainCount = SIZE(m_Chains);
	uint DBSize = ChainCount;
	if (optset_dbsize)
		DBSize = uint(opt_dbsize);
	ProgressLog("dbsize=%u\n", DBSize);

	DSSParams Params;
	Params.SetFromCmdLine(DBSize);

	DSS D;
	D.SetParams(Params);

	DSSAligner DA;
	DA.m_Params = &Params;

	m_Profiles.resize(ChainCount);
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "m_Profiles");
		D.Init(*m_Chains[i]);
		vector<vector<byte> > *Profile = new vector<vector<byte> >;
		D.GetProfile(*Profile);
		m_Profiles[i] = Profile;
		}

	m_HSMinScore = opt_hsminscore;
	m_MinDiagScore = opt_mindiagscore;
	const string &RefFN = opt_ref; // "c:/int/reseek_bench/alns/reseek_fast.tsv";

	m_fOut = CreateStdioFile(opt_output);

	Progress("Read Mu fasta...");
	SeqDB Input;
	Input.FromFasta(g_Arg1);
	m_Input = &Input;
	m_Input->SetLabelToIndex();
	Progress(" done.\n");
	m_SeqCount = m_Input->GetSeqCount();
	m_MuLettersVec.resize(m_SeqCount);

	m_MuSeqs = myalloc(byte *, m_SeqCount);
	for (uint SeqIdx = 0; SeqIdx < m_SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, m_SeqCount, "Mu letters");
		const string &Label = m_Input->GetLabel(SeqIdx);
		const string &ChainLabel = m_Chains[SeqIdx]->m_Label;
		asserta(Label == ChainLabel);
		uint L = m_Input->GetSeqLength(SeqIdx);
		byte *MuSeq = myalloc(byte, L);
		StrToMuLetters(m_Input->GetSeq(SeqIdx), MuSeq);
		m_MuSeqs[SeqIdx] = MuSeq;

		vector<byte> &MuLetters = m_MuLettersVec[SeqIdx];
		MuLetters.reserve(L);
		for (uint i = 0; i < L; ++i)
			MuLetters.push_back(MuSeq[i]);
		}

	m_MD.FromSeqDB(*m_Input);
	m_k = m_MD.m_k;
	m_DictSize = m_MD.m_DictSize;

	extern const short * const *Mu_S_ij_short;
	m_MM.Init(Mu_S_ij_short, 5, 36, 2);
	asserta(m_MM.m_AS_pow[5] == m_DictSize);

	DiagHSP DH;

// Buffer for high-scoring m_k-mers
	uint *HSKmers = myalloc(uint, m_DictSize);

	int *SeqIdxTToBestDiagScore = myalloc(int, m_SeqCount);
	int *SeqIdxTToBestDiag = myalloc(int, m_SeqCount);
	//int *SeqIdxTToBestDiagLo = myalloc(int, m_SeqCount);
	//int *SeqIdxTToBestDiagLen = myalloc(int, m_SeqCount);
	uint *SeqIdxTs = myalloc(uint, m_SeqCount);
	uint HitCount = 0;

	zero_array(SeqIdxTToBestDiagScore, m_SeqCount);
	zero_array(SeqIdxTToBestDiag, m_SeqCount);
	//zero_array(SeqIdxTToBestDiagLo, m_SeqCount);
	//zero_array(SeqIdxTToBestDiagLen, m_SeqCount);
	zero_array(SeqIdxTs, m_SeqCount);

	TwoHitDiag TH;
	for (uint SeqIdxQ = 0; SeqIdxQ < m_SeqCount; ++SeqIdxQ)
		{
		DoQ(SeqIdxQ, DA, TH, DH, HSKmers, SeqIdxTToBestDiagScore,
			SeqIdxTToBestDiag, SeqIdxTs);
		}
	CloseStdioFile(m_fOut);

	double PctIdeal = m_NT3*100.0/173058;

	ProgressLog("m_NT2=%u", m_NT2);
	ProgressLog("\tNT3=%u", m_NT3);
	ProgressLog("\tNH2=%u", m_NH2);
	ProgressLog("\tIdeal=%.1f%%", PctIdeal);
	ProgressLog("\ths=%u", opt_hsminscore);
	ProgressLog("\tdg=%u", opt_mindiagscore);
	}

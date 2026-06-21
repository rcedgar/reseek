#include "myutils.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_chain_reader.h"
#include "bcadata.h"
#include "chaq.h"
#include <thread>
#include <set>

static FILE *s_fCal = 0;
static FILE *s_fCan = 0;
static FILE *s_fFasta = 0;
static FILE *s_fNuHexFasta = 0;
static FILE *s_fKappaFasta = 0;

static PDBFileScanner *s_ptrFS = 0;
static BCAData *s_ptrBCA = 0;
static BCAData *s_ptrBCB = 0;

static bool s_want_bcb = false;
static bool s_want_kappafasta = false;

static mutex s_LockStats;
static mutex s_LockCal;
static mutex s_LockCan;
static mutex s_LockFasta;
static mutex s_LockNuHexFasta;
static mutex s_LockKappaFasta;
static mutex s_LockBCA;

static uint s_MinChainLength = 1;
static uint s_Converted = 0;
static uint s_TooShort = 0;
static uint s_Shortest = UINT_MAX;
static uint s_InputCount = 0;
static uint s_OutputCount = 0;
static time_t s_LastTime = 0;
static set<string> *s_ptrLabelSet = 0;
static uint s_LabelSetSize = 0;

static void WriteCan(FILE *f, const flat_chain_t *chain)
	{
	if (f == 0)
		return;
	const uint L = chain->get_length();
	asserta(chain->has_nu());
	const uint8_t *nu = chain->get_nu_data();
	fprintf(f, ">%s\n", chain->m_label.c_str());
	for (uint pos = 0; pos < L; ++pos)
		{
		char aa = chain->get_aa(pos);
		float x, y, z;
		chain->get_coords(pos, x, y, z);
		fprintf(f, "%c\t%.1f\t%.1f\t%.1f\t%02x\n",
			aa, x, y, z, nu[pos]);
		}
	}

static void WriteKappaFasta(FILE *f, const flat_chain_t *chain,
	uint8_t *codeseq_kappa)
	{
	if (f == 0)
		return;
	const uint L = chain->get_length();
	asserta(chain->has_nu());
	const uint8_t *nu = chain->get_nu_data();
	chaq::codeseq_nu_to_kappa(
		nu, L, codeseq_kappa, flat_params::m_maxL);
	codeseq_to_fasta(f, chain->m_label, codeseq_kappa, L, KAPPA_AS);
	}

static void ThreadBody(uint ThreadIndex)
	{
	chaq_vecs2 cv;
	uint8_t *codeseq_kappa = 0;
	if (s_want_bcb)
		chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
	if (s_want_kappafasta)
		codeseq_kappa = myalloc(uint8_t, flat_params::m_maxL);

	flat_chain_reader CR;
	CR.Open(*s_ptrFS);

	for (;;)
		{
		s_LockStats.lock();
		time_t Now = time(0);
		if (Now - s_LastTime > 0)
			{
			if (s_LabelSetSize > 0)
				Progress("%u / %u chains found (%s searched)\r",
					s_OutputCount, s_LabelSetSize,
					IntToStr(s_InputCount));
			else if (s_TooShort > 0)
				Progress("%s chains, %.1f%% too short (min %u, shortest %u)\r",
				  IntToStr(s_Converted), GetPct(s_TooShort, s_Converted),
				  s_MinChainLength, s_Shortest);
			else
				Progress("%s chains converted\r", IntToStr(s_Converted));
			uint ne = flat_chain_reader::m_CRGlobalFormatErrors;
			if (ne > 0)
				Progress(" %u format errors", ne);
			s_LastTime = Now;
			}
		s_LockStats.unlock();

		flat_chain_t *chain = CR.GetNext();
		if (chain == 0)
			break;

		const uint L = chain->get_length();
		asserta(L > 0);
		asserta(L < flat_params::m_maxL);

		s_LockStats.lock();
		++s_InputCount;
		s_Shortest = min(L, s_Shortest);
		s_LockStats.unlock();

		if (s_ptrLabelSet != 0)
			{
			string UpperLabel = chain->m_label;
			ToUpper(UpperLabel);
			set<string>::iterator iter = s_ptrLabelSet->find(UpperLabel);
			if (iter == s_ptrLabelSet->end())
				{
				delete chain;
				continue;
				}
			}

		if (L < s_MinChainLength)
			{
			s_LockStats.lock();
			++s_TooShort;
			s_LockStats.unlock();
			delete chain;
			continue;
			}

		if (optset_subsample)
			{
			s_LockStats.lock();
			uint n = s_InputCount;
			s_LockStats.unlock();
			if (n%opt(subsample) != 0)
				{
				delete chain;
				continue;
				}
			}

		s_LockStats.lock();
		++s_OutputCount;
		if (s_ptrLabelSet != 0)
			{
			string UpperLabel = chain->m_label;
			ToUpper(UpperLabel);
			s_ptrLabelSet->erase(UpperLabel);
			}
		s_LockStats.unlock();

		asserta(chain->has_nu());
		const uint8_t *nu = chain->get_nu_data();

		if (s_fCal != 0)
			{
			s_LockCal.lock();
			chain->to_cal(s_fCal);
			s_LockCal.unlock();
			}

		if (s_fCan != 0)
			{
			s_LockCan.lock();
			WriteCan(s_fCan, chain);
			s_LockCan.unlock();
			}

		if (s_fFasta != 0)
			{
			s_LockFasta.lock();
			chain->to_fasta(s_fFasta);
			s_LockFasta.unlock();
			}

		if (s_fNuHexFasta != 0)
			{
			s_LockNuHexFasta.lock();
			codeseq_to_hexfasta(s_fNuHexFasta, chain->m_label, nu, L);
			s_LockNuHexFasta.unlock();
			}

		if (s_fKappaFasta != 0)
			{
			s_LockKappaFasta.lock();
			WriteKappaFasta(s_fKappaFasta, chain, codeseq_kappa);
			s_LockKappaFasta.unlock();
			}

		if (s_ptrBCA != 0 || s_ptrBCB != 0)
			{
			s_LockBCA.lock();
			if (s_ptrBCA != 0)
				s_ptrBCA->write_flat_chain(chain, &cv);
			if (s_ptrBCB != 0)
				s_ptrBCB->write_flat_chain(chain, &cv);
			s_LockBCA.unlock();
			}

		s_LockStats.lock();
		++s_Converted;
		s_LockStats.unlock();

		delete chain;
		}

	if (s_want_bcb)
		chaq::free_chaq_vecs2(cv);
	myfree(codeseq_kappa);
	}

void cmd_flat_convert()
	{
	if (optset_output)
		Die("Use -cal, -can, -bca, -bcb, -fasta, -nuhexfasta or "
		  "-kappafasta not -output");

	const bool want_cal = optset_cal;
	const bool want_can = optset_can;
	const bool want_bca = optset_bca;
	s_want_bcb = optset_bcb;
	const bool want_fasta = optset_fasta;
	const bool want_nuhexfasta = optset_nuhexfasta;
	s_want_kappafasta = optset_kappafasta;

	if (!want_cal && !want_can && !want_bca && !s_want_bcb &&
		!want_fasta && !want_nuhexfasta && !s_want_kappafasta)
		Die("Must specify one or more output options: "
		  "-cal, -can, -bca, -bcb, -fasta, -nuhexfasta, -kappafasta");

	s_MinChainLength = 1;
	if (optset_minchainlength)
		s_MinChainLength = opt(minchainlength);

	vector<string> Labels;
	if (optset_labels)
		{
		ReadLinesFromFile(opt(labels), Labels);
		if (Labels.empty())
			Die("No labels found in '%s'", opt(labels));
		}
	else if (optset_label)
		Labels.push_back(opt(label));

	set<string> LabelSet;
	for (uint i = 0; i < SIZE(Labels); ++i)
		{
		ToUpper(Labels[i]);
		LabelSet.insert(Labels[i]);
		}
	if (!LabelSet.empty())
		{
		s_ptrLabelSet = &LabelSet;
		s_LabelSetSize = SIZE(LabelSet);
		}
	else
		{
		s_ptrLabelSet = 0;
		s_LabelSetSize = 0;
		}

	PDBFileScanner FS;
	FS.Open(g_Arg1);
	s_ptrFS = &FS;

	s_fCal = CreateStdioFile(opt(cal));
	s_fCan = CreateStdioFile(opt(can));
	s_fFasta = CreateStdioFile(opt(fasta));
	s_fNuHexFasta = CreateStdioFile(opt(nuhexfasta));
	s_fKappaFasta = CreateStdioFile(opt(kappafasta));

	BCAData BCA;
	BCAData BCB;
	s_ptrBCA = 0;
	s_ptrBCB = 0;
	if (want_bca)
		{
		BCA.Create(opt(bca), false);
		s_ptrBCA = &BCA;
		}
	if (s_want_bcb)
		{
		BCB.Create(opt(bcb), true);
		s_ptrBCB = &BCB;
		}

	s_InputCount = 0;
	s_Converted = 0;
	s_OutputCount = 0;
	s_TooShort = 0;
	s_Shortest = UINT_MAX;
	s_LastTime = 0;

	vector<thread *> ts;
	const uint ThreadCount = GetRequestedThreadCount();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts.push_back(new thread(ThreadBody, ThreadIndex));
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		ts[ThreadIndex]->join();
		delete ts[ThreadIndex];
		}

	if (s_InputCount > 10000)
		{
		ProgressLog("\n");
		ProgressLogPrefix("%u / %u converted (%s, %.1f%%)\n",
		  s_OutputCount, s_InputCount, IntToStr(s_InputCount),
		  GetPct(s_OutputCount, s_InputCount));
		if (s_TooShort > 0)
			ProgressLogPrefix("%u too short (%s, %.1f%%) min length %u\n",
			  s_TooShort, IntToStr(s_TooShort),
			  GetPct(s_TooShort, s_InputCount), s_MinChainLength);
		}
	else
		{
		ProgressLog("\n");
		ProgressLogPrefix("%u converted\n", s_InputCount);
		if (s_TooShort > 0)
			ProgressLogPrefix("%u too short (%.1f%%), min length %u, shortest %u\n",
			  s_TooShort, GetPct(s_TooShort, s_InputCount),
			  s_MinChainLength, s_Shortest);
		}
	uint ne = flat_chain_reader::m_CRGlobalFormatErrors;
	if (ne > 0)
		ProgressLogPrefix("%u format errors\n", ne);

	if (optset_label || optset_labels)
		{
		ProgressLogPrefix("Searched for %u labels, %u found (%.1f%%)\n",
			s_LabelSetSize, s_OutputCount,
			GetPct(s_OutputCount, s_LabelSetSize));
		uint NotFound = SIZE(*s_ptrLabelSet);
		if (NotFound > 0)
			{
			Progress("%u not found", NotFound);
			Log("%u not found\n", NotFound);
			uint Counter = 0;
			for (set<string>::const_iterator iter = s_ptrLabelSet->begin();
				 iter != s_ptrLabelSet->end(); ++iter)
				{
				if (Counter++ < 3)
					Progress(" %s", iter->c_str());
				Log(">%s\n", iter->c_str());
				}
			Progress("\n");
			}
		}

	CloseStdioFile(s_fCal);
	CloseStdioFile(s_fCan);
	CloseStdioFile(s_fFasta);
	CloseStdioFile(s_fNuHexFasta);
	CloseStdioFile(s_fKappaFasta);
	s_fCal = 0;
	s_fCan = 0;
	s_fFasta = 0;
	s_fNuHexFasta = 0;
	s_fKappaFasta = 0;

	if (want_bca)
		{
		Progress("finalizing BCA... ");
		BCA.Close();
		Progress("done\n");
		}
	if (s_want_bcb)
		{
		Progress("finalizing BCB... ");
		BCB.Close();
		Progress("done\n");
		}

	s_ptrBCA = 0;
	s_ptrBCB = 0;
	s_ptrFS = 0;
	s_ptrLabelSet = 0;
	}

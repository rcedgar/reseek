#include "myutils.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_chain_reader.h"
#include "bcadata.h"
#include "chaq.h"
#include <thread>
#include <atomic>
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
static bool s_want_pdbcaoutdir = false;
static string s_pdbcaoutdir;

static mutex s_LockStats;
static mutex s_LockCal;
static mutex s_LockCan;
static mutex s_LockFasta;
static mutex s_LockNuHexFasta;
static mutex s_LockKappaFasta;
static mutex s_LockBCA;
static mutex s_LockPdbCaOut;

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

static string MakeUniquePdbPath(const string &dir, const string &label)
	{
	string base = dir;
	Dirize(base);
	base += label;
	string path = base + ".pdb";
	if (!StdioFileExists(path))
		return path;
	for (uint n = 1; n < 10000; ++n)
		{
		string pathN;
		Ps(pathN, "%s_dupe%03u.pdb", base.c_str(), n);
		if (!StdioFileExists(pathN))
			return pathN;
		}
	Die("Too many duplicate PDB files for label '%s'", label.c_str());
	return "";
	}

static void WritePdbCa(const flat_chain_t *chain)
	{
	s_LockPdbCaOut.lock();
	string label = chain->m_label;
	trunc_label(label);
	string path = MakeUniquePdbPath(s_pdbcaoutdir, label);
	char chainId = ExtractChainIdFromLabel(label);
	chain->to_pdb(path, chainId);
	s_LockPdbCaOut.unlock();
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
		asserta(L <= flat_params::m_maxL);

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

		if (s_want_pdbcaoutdir)
			WritePdbCa(chain);

		s_LockStats.lock();
		++s_Converted;
		s_LockStats.unlock();

		delete chain;
		}

	if (s_want_bcb)
		chaq::free_chaq_vecs2(cv);
	myfree(codeseq_kappa);
	}

// ---- Fast parallel path: single .bca/.bcb source -> -bca/-bcb output -----
//
// The generic path routes a single input file through one PDBFileScanner,
// so only one thread does any work. For an indexed source (.bca/.bcb) we can
// instead partition chains by index and have every thread read + compute nu
// + write its own output shard in parallel. Limited to BCA/BCB outputs with
// no label filtering; everything else falls back to the generic path.

static BCAData *s_fast_src = 0;
static BCAData *s_fast_bca = 0;
static BCAData *s_fast_bcb = 0;
static uint s_fast_N = 0;
static uint s_fast_minlen = 1;
static std::atomic<uint> s_fast_next;
static std::atomic<uint> s_fast_done;
static mutex s_fast_stats_lock;
static uint s_fast_input = 0;
static uint s_fast_converted = 0;
static uint s_fast_tooshort = 0;
static uint s_fast_shortest = UINT_MAX;

static void FastThreadBody(uint ThreadIndex)
	{
	const uint maxL = flat_params::m_maxL;
	uint8_t *nu_read = 0;
	if (s_fast_src->m_HasNuSequences)
		nu_read = myalloc(uint8_t, maxL);

	uint loc_input = 0;
	uint loc_conv = 0;
	uint loc_short = 0;
	uint loc_shortest = UINT_MAX;

	const uint CHUNK = 256;
	time_t LastTime = 0;
	for (;;)
		{
		uint begin = s_fast_next.fetch_add(CHUNK);
		if (begin >= s_fast_N)
			break;
		uint end = begin + CHUNK;
		if (end > s_fast_N)
			end = s_fast_N;

		for (uint idx = begin; idx < end; ++idx)
			{
			flat_chain_t *chain = s_fast_src->read_flat_chain(idx);
			const uint L = chain->get_length();
			if (L == 0)
				{
				delete chain;
				continue;
				}
			++loc_input;
			loc_shortest = min(L, loc_shortest);

			if (L < s_fast_minlen)
				{
				++loc_short;
				delete chain;
				continue;
				}

			if (optset_subsample && ((idx + 1) % opt(subsample)) != 0)
				{
				delete chain;
				continue;
				}

// Carry forward stored nu for BCB sources (read_codeseq_nu is not internally
// locked, so guard the shared source handle here).
			if (nu_read != 0)
				{
				s_fast_src->m_ReadLock.lock();
				uint nL = s_fast_src->read_codeseq_nu(nu_read, idx, maxL);
				s_fast_src->m_ReadLock.unlock();
				asserta(nL == L);
				chain->set_nu_codes(nu_read, L);
				}

			if (s_fast_bca != 0)
				s_fast_bca->write_flat_chain_shard(ThreadIndex, chain);
			if (s_fast_bcb != 0)
				s_fast_bcb->write_flat_chain_shard(ThreadIndex, chain);
			++loc_conv;
			delete chain;
			}

		s_fast_done.fetch_add(end - begin);
		if (ThreadIndex == 0)
			{
			time_t Now = time(0);
			if (Now - LastTime > 0)
				{
				uint done = s_fast_done.load();
				Progress("%u / %u chains converted (%.2f%%)\r",
					done, s_fast_N, GetPct(done, s_fast_N));
				LastTime = Now;
				}
			}
		}

	myfree(nu_read);

	s_fast_stats_lock.lock();
	s_fast_input += loc_input;
	s_fast_converted += loc_conv;
	s_fast_tooshort += loc_short;
	if (loc_shortest < s_fast_shortest)
		s_fast_shortest = loc_shortest;
	s_fast_stats_lock.unlock();
	}

static bool FastPathEligible(bool want_cal, bool want_can, bool want_bca,
	bool want_bcb, bool want_fasta, bool want_nuhexfasta,
	bool want_kappafasta, bool want_pdbcaoutdir, bool have_labels)
	{
	if (want_cal || want_can || want_fasta || want_nuhexfasta ||
		want_kappafasta || want_pdbcaoutdir)
		return false;
	if (!(want_bca || want_bcb))
		return false;
	if (have_labels)
		return false;
	if (!IsRegularFile(g_Arg1))
		return false;
	string Ext;
	GetExtFromPathName(g_Arg1, Ext);
	ToLower(Ext);
	return (Ext == "bca" || Ext == "bcb");
	}

static void RunFastBcx(bool want_bca, bool want_bcb)
	{
	BCAData src;
	src.Open(g_Arg1);
	const uint N = src.GetChainCount();

	BCAData out_bca;
	BCAData out_bcb;
	const uint ThreadCount = GetRequestedThreadCount();
	if (want_bca)
		out_bca.CreateSharded(opt(bca), false, ThreadCount);
	if (want_bcb)
		out_bcb.CreateSharded(opt(bcb), true, ThreadCount);

	s_fast_src = &src;
	s_fast_bca = (want_bca ? &out_bca : 0);
	s_fast_bcb = (want_bcb ? &out_bcb : 0);
	s_fast_N = N;
	s_fast_minlen = s_MinChainLength;
	s_fast_next = 0;
	s_fast_done = 0;
	s_fast_input = 0;
	s_fast_converted = 0;
	s_fast_tooshort = 0;
	s_fast_shortest = UINT_MAX;

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts.push_back(new thread(FastThreadBody, ThreadIndex));
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		ts[ThreadIndex]->join();
		delete ts[ThreadIndex];
		}

	if (want_bca)
		{
		Progress("finalizing BCA... ");
		out_bca.CloseSharded();
		Progress("done\n");
		}
	if (want_bcb)
		{
		Progress("finalizing BCB... ");
		out_bcb.CloseSharded();
		Progress("done\n");
		}
	src.Close();

	ProgressLog("\n");
	ProgressLogPrefix("%u / %u converted (%s)\n",
		s_fast_converted, s_fast_input, IntToStr(s_fast_input));
	if (s_fast_tooshort > 0)
		ProgressLogPrefix("%u too short (min length %u, shortest %u)\n",
			s_fast_tooshort, s_fast_minlen,
			(s_fast_shortest == UINT_MAX ? 0 : s_fast_shortest));

	s_fast_src = 0;
	s_fast_bca = 0;
	s_fast_bcb = 0;
	}

void cmd_flat_convert()
	{
	if (optset_output)
		Die("Use -cal, -can, -bca, -bcb, -fasta, -nuhexfasta, "
		  "-kappafasta or -pdbcaoutdir not -output");

	const bool want_cal = optset_cal;
	const bool want_can = optset_can;
	const bool want_bca = optset_bca;
	s_want_bcb = optset_bcb;
	const bool want_fasta = optset_fasta;
	const bool want_nuhexfasta = optset_nuhexfasta;
	s_want_kappafasta = optset_kappafasta;
	s_want_pdbcaoutdir = optset_pdbcaoutdir;
	s_pdbcaoutdir = opt(pdbcaoutdir);

	if (!want_cal && !want_can && !want_bca && !s_want_bcb &&
		!want_fasta && !want_nuhexfasta && !s_want_kappafasta &&
		!s_want_pdbcaoutdir)
		Die("Must specify one or more output options: "
		  "-cal, -can, -bca, -bcb, -fasta, -nuhexfasta, -kappafasta, "
		  "-pdbcaoutdir");

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

	if (FastPathEligible(want_cal, want_can, want_bca, s_want_bcb,
		want_fasta, want_nuhexfasta, s_want_kappafasta,
		s_want_pdbcaoutdir, s_ptrLabelSet != 0))
		{
		RunFastBcx(want_bca, s_want_bcb);
		uint ne = flat_chain_reader::m_CRGlobalFormatErrors;
		if (ne > 0)
			ProgressLogPrefix("%u format errors\n", ne);
		s_ptrFS = 0;
		s_ptrLabelSet = 0;
		return;
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

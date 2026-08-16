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

// ---- Fast parallel path: single .bca/.bcb/.cal source -> text and/or bcx -----
//
// The generic path routes a single input file through one PDBFileScanner,
// so only one thread does any work. For an indexed source (.bca/.bcb) we
// partition chains by index. For a text .cal we split the file by byte
// range; each thread skips to the next record start and parses in place.
// Threads write their own output shard(s). Text shards are concatenated
// in thread order at close. No label filtering; -pdbcaoutdir stays on
// the generic path.

struct text_shards_t
	{
	string final_fn;
	vector<string> tmp_fns;
	vector<FILE *> fs;

	void open(const string &fn, uint n)
		{
		asserta(fn != "");
		asserta(n >= 1);
		final_fn = fn;
		tmp_fns.resize(n);
		fs.resize(n);
		for (uint i = 0; i < n; ++i)
			{
			Ps(tmp_fns[i], "%s.tmp.%u", fn.c_str(), i);
			fs[i] = CreateStdioFile(tmp_fns[i]);
			}
		}

	FILE *get(uint tid) const
		{
		asserta(tid < SIZE(fs));
		return fs[tid];
		}

	void close_files()
		{
		for (uint i = 0; i < SIZE(fs); ++i)
			CloseStdioFile(fs[i]);
		fs.clear();
		}
	};

static void ConcatShardFiles(const string &dest, const vector<string> &shards)
	{
	FILE *fout = CreateStdioFile(dest);
	const uint BUF = 1u << 20;
	byte *buf = myalloc(byte, BUF);
	for (uint i = 0; i < SIZE(shards); ++i)
		{
		FILE *fin = OpenStdioFile(shards[i]);
		uint64 sz = GetStdioFileSize64(fin);
		uint64 pos = 0;
		while (pos < sz)
			{
			uint64 n = sz - pos;
			if (n > BUF)
				n = BUF;
			ReadStdioFile64(fin, pos, buf, n);
			WriteStdioFile64(fout, buf, n);
			pos += n;
			}
		CloseStdioFile(fin);
		DeleteStdioFile(shards[i]);
		}
	CloseStdioFile(fout);
	myfree(buf);
	}

static void FinalizeTextShards(text_shards_t &sh, const char *what)
	{
	if (sh.final_fn == "")
		return;
	sh.close_files();
	Progress("finalizing %s... ", what);
	ConcatShardFiles(sh.final_fn, sh.tmp_fns);
	Progress("done\n");
	sh.tmp_fns.clear();
	sh.final_fn.clear();
	}

static BCAData *s_fast_src = 0;
static BCAData *s_fast_bca = 0;
static BCAData *s_fast_bcb = 0;
static text_shards_t *s_fast_cal = 0;
static text_shards_t *s_fast_can = 0;
static text_shards_t *s_fast_fasta = 0;
static text_shards_t *s_fast_nuhex = 0;
static text_shards_t *s_fast_kappa = 0;
static bool s_fast_need_nu = false;
static bool s_fast_need_nu_on_chain = false;
static uint s_fast_N = 0;
static uint s_fast_minlen = 1;
static std::atomic<uint> s_fast_next;
static std::atomic<uint> s_fast_done;
static mutex s_fast_stats_lock;
static uint s_fast_input = 0;
static uint s_fast_converted = 0;
static uint s_fast_tooshort = 0;
static uint s_fast_shortest = UINT_MAX;

static FILE *s_fast_cal_in = 0;
static uint64 s_fast_cal_insz = 0;
static uint s_fast_cal_nthreads = 1;
static const uint CAL_FAST_BUF = 4u << 20;

static void FastThreadBody(uint ThreadIndex)
	{
	const uint maxL = flat_params::m_maxL;
	uint8_t *nu_buf = 0;
	uint8_t *codeseq_kappa = 0;
	sid_t *distmx = 0;
	chaq_vecs2 cv;
	bool cv_inited = false;

	if (s_fast_need_nu)
		nu_buf = myalloc(uint8_t, maxL);
	if (s_fast_kappa != 0)
		codeseq_kappa = myalloc(uint8_t, maxL);

	FILE *fCal = (s_fast_cal ? s_fast_cal->get(ThreadIndex) : 0);
	FILE *fCan = (s_fast_can ? s_fast_can->get(ThreadIndex) : 0);
	FILE *fFasta = (s_fast_fasta ? s_fast_fasta->get(ThreadIndex) : 0);
	FILE *fNuHex = (s_fast_nuhex ? s_fast_nuhex->get(ThreadIndex) : 0);
	FILE *fKappa = (s_fast_kappa ? s_fast_kappa->get(ThreadIndex) : 0);

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

			if (s_fast_need_nu_on_chain)
				{
				if (s_fast_src->m_HasNuSequences)
					{
					uint nL = s_fast_src->read_codeseq_nu(nu_buf, idx, maxL);
					asserta(nL == L);
					chain->set_nu_codes(nu_buf, L);
					}
				else
					{
					if (!cv_inited)
						{
						distmx = myalloc(sid_t,
							flat_params::m_distmx_bandwidth*maxL);
						chaq::alloc_chaq_vecs2(cv, maxL);
						cv_inited = true;
						}
					chaq::fill_codeseq_nu_from_chain(
						chain, distmx, &cv, nu_buf, maxL);
					chain->set_nu_codes(nu_buf, L);
					}
				}
			else if (s_fast_need_nu && s_fast_src->m_HasNuSequences &&
				s_fast_bcb != 0)
				{
// Pass stored nu through to BCB shard writer (avoids recomputing).
				uint nL = s_fast_src->read_codeseq_nu(nu_buf, idx, maxL);
				asserta(nL == L);
				chain->set_nu_codes(nu_buf, L);
				}

			if (fCal != 0)
				chain->to_cal(fCal);
			if (fCan != 0)
				WriteCan(fCan, chain);
			if (fFasta != 0)
				chain->to_fasta(fFasta);
			if (fNuHex != 0)
				{
				asserta(chain->has_nu());
				codeseq_to_hexfasta(fNuHex, chain->m_label,
					chain->get_nu_data(), L);
				}
			if (fKappa != 0)
				WriteKappaFasta(fKappa, chain, codeseq_kappa);
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

	myfree(nu_buf);
	myfree(codeseq_kappa);
	if (cv_inited)
		{
		myfree(distmx);
		chaq::free_chaq_vecs2(cv);
		}

	s_fast_stats_lock.lock();
	s_fast_input += loc_input;
	s_fast_converted += loc_conv;
	s_fast_tooshort += loc_short;
	if (loc_shortest < s_fast_shortest)
		s_fast_shortest = loc_shortest;
	s_fast_stats_lock.unlock();
	}

struct cal_buf_t
	{
	FILE *f = 0;
	char *buf = 0;
	uint n = 0;
	uint i = 0;
	uint64 buf_pos = 0;
	uint64 file_size = 0;
	const char *fn = 0;
	bool pending = false;
	uint64 pending_pos = 0;
	string pending_label;

	uint64 tell() const { return buf_pos + i; }

	void open(FILE *F, uint64 start, uint64 fsz, const char *FN)
		{
		f = F;
		file_size = fsz;
		fn = FN;
		buf = myalloc(char, CAL_FAST_BUF);
		buf_pos = start;
		i = 0;
		n = 0;
		refill();
		}

	void close_buf()
		{
		myfree(buf);
		buf = 0;
		}

	void refill()
		{
		if (i > 0)
			{
			uint rem = n - i;
			if (rem > 0)
				memmove(buf, buf + i, rem);
			buf_pos += i;
			n = rem;
			i = 0;
			}
		if (n == CAL_FAST_BUF)
			return;
		uint64 pos = buf_pos + n;
		if (pos >= file_size)
			return;
		uint64 room = CAL_FAST_BUF - n;
		uint64 left = file_size - pos;
		if (room > left)
			room = left;
		n += uint(ReadStdioFile64_NoFail(f, pos, buf + n, room));
		}

	bool at_eof()
		{
		if (i < n)
			return false;
		refill();
		return i >= n;
		}

	bool read_line(const char *&p, uint &len)
		{
		if (at_eof())
			return false;
		uint start = i;
		for (;;)
			{
			while (i < n && buf[i] != '\n')
				++i;
			if (i < n)
				break;
			if (start == 0 && n == CAL_FAST_BUF)
				Die("%s: CAL line exceeds %u bytes", fn, CAL_FAST_BUF);
			i = start;
			uint have = n - i;
			refill();
			start = 0;
			i = have;
			if (i >= n)
				{
				p = buf + start;
				len = n - start;
				if (len > 0 && p[len - 1] == '\r')
					--len;
				i = n;
				return true;
				}
			}
		p = buf + start;
		len = i - start;
		if (len > 0 && p[len - 1] == '\r')
			--len;
		++i;
		return true;
		}
	};

static bool ParseCalResidueLine(const char *p, uint len,
	char &aa, float &x, float &y, float &z)
	{
	if (len < 7)
		return false;
	aa = p[0];
	if (p[1] != '\t')
		return false;
	const char *end = p + len;
	char *e = 0;
	const char *s = p + 2;
	x = strtof(s, &e);
	if (e == s || e >= end || *e != '\t')
		return false;
	s = e + 1;
	y = strtof(s, &e);
	if (e == s || e >= end || *e != '\t')
		return false;
	s = e + 1;
	z = strtof(s, &e);
	if (e == s)
		return false;
	while (e < end && (*e == ' ' || *e == '\t'))
		++e;
	return e == end;
	}

static void EmitFastCalChain(uint ThreadIndex, flat_chain_t *chain, uint idx,
	uint &loc_input, uint &loc_conv, uint &loc_short, uint &loc_shortest,
	FILE *fCal, FILE *fCan, FILE *fFasta, FILE *fNuHex, FILE *fKappa,
	uint8_t *nu_buf, uint8_t *codeseq_kappa,
	sid_t *&distmx, chaq_vecs2 &cv, bool &cv_inited)
	{
	const uint maxL = flat_params::m_maxL;
	const uint L = chain->get_length();
	if (L == 0)
		{
		delete chain;
		return;
		}
	++loc_input;
	loc_shortest = min(L, loc_shortest);

	if (L < s_fast_minlen)
		{
		++loc_short;
		delete chain;
		return;
		}

	if (optset_subsample && ((idx + 1) % opt(subsample)) != 0)
		{
		delete chain;
		return;
		}

	if (s_fast_need_nu_on_chain)
		{
		if (!cv_inited)
			{
			distmx = myalloc(sid_t,
				flat_params::m_distmx_bandwidth*maxL);
			chaq::alloc_chaq_vecs2(cv, maxL);
			cv_inited = true;
			}
		chaq::fill_codeseq_nu_from_chain(
			chain, distmx, &cv, nu_buf, maxL);
		chain->set_nu_codes(nu_buf, L);
		}

	if (fCal != 0)
		chain->to_cal(fCal);
	if (fCan != 0)
		WriteCan(fCan, chain);
	if (fFasta != 0)
		chain->to_fasta(fFasta);
	if (fNuHex != 0)
		{
		asserta(chain->has_nu());
		codeseq_to_hexfasta(fNuHex, chain->m_label,
			chain->get_nu_data(), L);
		}
	if (fKappa != 0)
		WriteKappaFasta(fKappa, chain, codeseq_kappa);
	if (s_fast_bca != 0)
		s_fast_bca->write_flat_chain_shard(ThreadIndex, chain);
	if (s_fast_bcb != 0)
		s_fast_bcb->write_flat_chain_shard(ThreadIndex, chain);
	++loc_conv;
	delete chain;
	}

static flat_chain_t *ParseNextCalChain(cal_buf_t &R, uint64 range_end,
	bool &skipping)
	{
	string label;
	bool in_chain = false;

	if (R.pending)
		{
		if (R.pending_pos >= range_end)
			return 0;
		label.swap(R.pending_label);
		R.pending = false;
		in_chain = true;
		skipping = false;
		}

	vector<char> aas;
	vector<float> Xs, Ys, Zs;
	aas.reserve(RESERVE_CHAIN_LENGTH);
	Xs.reserve(RESERVE_CHAIN_LENGTH);
	Ys.reserve(RESERVE_CHAIN_LENGTH);
	Zs.reserve(RESERVE_CHAIN_LENGTH);

	const uint maxL = flat_params::m_maxL;
	bool truncated = false;

	for (;;)
		{
		uint64 pos = R.tell();
		const char *p = 0;
		uint len = 0;
		if (!R.read_line(p, len))
			break;
		if (len == 0)
			continue;

		if (p[0] == '>')
			{
			if (in_chain)
				{
				R.pending = true;
				R.pending_pos = pos;
				R.pending_label.assign(p + 1, len - 1);
				break;
				}
			if (pos >= range_end)
				return 0;
			if (skipping)
				skipping = false;
			label.assign(p + 1, len - 1);
			in_chain = true;
			continue;
			}

		if (!in_chain)
			{
			if (skipping)
				continue;
			Die("%s: Expected '>' in CAL file", R.fn);
			}

		if (SIZE(aas) >= maxL)
			{
			truncated = true;
			continue;
			}

		char aa;
		float x, y, z;
		if (!ParseCalResidueLine(p, len, aa, x, y, z))
			Die("%s: Invalid CAL record '%.*s'", R.fn, int(len), p);

		aas.push_back(aa);
		Xs.push_back(x);
		Ys.push_back(y);
		Zs.push_back(z);
		}

	if (!in_chain)
		return 0;
	if (truncated)
		g_flat_n_truncated_chains.fetch_add(1, memory_order_relaxed);
	if (aas.empty())
		return flat_chain_t::newflat(0);

	return flat_chain_t::newflat(label, aas, Xs, Ys, Zs);
	}

static void FastCalThreadBody(uint ThreadIndex)
	{
	const uint maxL = flat_params::m_maxL;
	uint8_t *nu_buf = 0;
	uint8_t *codeseq_kappa = 0;
	sid_t *distmx = 0;
	chaq_vecs2 cv;
	bool cv_inited = false;

	if (s_fast_need_nu_on_chain)
		nu_buf = myalloc(uint8_t, maxL);
	if (s_fast_kappa != 0)
		codeseq_kappa = myalloc(uint8_t, maxL);

	FILE *fCal = (s_fast_cal ? s_fast_cal->get(ThreadIndex) : 0);
	FILE *fCan = (s_fast_can ? s_fast_can->get(ThreadIndex) : 0);
	FILE *fFasta = (s_fast_fasta ? s_fast_fasta->get(ThreadIndex) : 0);
	FILE *fNuHex = (s_fast_nuhex ? s_fast_nuhex->get(ThreadIndex) : 0);
	FILE *fKappa = (s_fast_kappa ? s_fast_kappa->get(ThreadIndex) : 0);

	uint loc_input = 0;
	uint loc_conv = 0;
	uint loc_short = 0;
	uint loc_shortest = UINT_MAX;

	const uint T = s_fast_cal_nthreads;
	const uint64 sz = s_fast_cal_insz;
	const uint64 range_begin = (sz * ThreadIndex) / T;
	const uint64 range_end = (sz * (ThreadIndex + 1)) / T;

	cal_buf_t R;
	R.open(s_fast_cal_in, range_begin, sz, g_Arg1.c_str());
	bool skipping = (range_begin > 0);

	time_t LastTime = 0;
	for (;;)
		{
		flat_chain_t *chain = ParseNextCalChain(R, range_end, skipping);
		if (chain == 0)
			break;

		uint idx = s_fast_next.fetch_add(1);
		EmitFastCalChain(ThreadIndex, chain, idx,
			loc_input, loc_conv, loc_short, loc_shortest,
			fCal, fCan, fFasta, fNuHex, fKappa,
			nu_buf, codeseq_kappa, distmx, cv, cv_inited);

		uint done = s_fast_done.fetch_add(1) + 1;
		if (ThreadIndex == 0)
			{
			time_t Now = time(0);
			if (Now - LastTime > 0)
				{
				Progress("%s chains converted\r", FloatToStr(done));
				LastTime = Now;
				}
			}
		}

	R.close_buf();
	myfree(nu_buf);
	myfree(codeseq_kappa);
	if (cv_inited)
		{
		myfree(distmx);
		chaq::free_chaq_vecs2(cv);
		}

	s_fast_stats_lock.lock();
	s_fast_input += loc_input;
	s_fast_converted += loc_conv;
	s_fast_tooshort += loc_short;
	if (loc_shortest < s_fast_shortest)
		s_fast_shortest = loc_shortest;
	s_fast_stats_lock.unlock();
	}

static void RunFastCal(bool want_cal, bool want_can, bool want_bca,
	bool want_bcb, bool want_fasta, bool want_nuhexfasta,
	bool want_kappafasta)
	{
	FILE *fin = OpenStdioFile(g_Arg1);
	const uint64 sz = GetStdioFileSize64(fin);
	if (sz == 0)
		Die("Empty CAL file '%s'", g_Arg1.c_str());

	BCAData out_bca;
	BCAData out_bcb;
	text_shards_t out_cal;
	text_shards_t out_can;
	text_shards_t out_fasta;
	text_shards_t out_nuhex;
	text_shards_t out_kappa;
	const uint ThreadCount = GetRequestedThreadCount();
	if (want_bca)
		out_bca.CreateSharded(opt(bca), false, ThreadCount);
	if (want_bcb)
		out_bcb.CreateSharded(opt(bcb), true, ThreadCount);
	if (want_cal)
		out_cal.open(opt(cal), ThreadCount);
	if (want_can)
		out_can.open(opt(can), ThreadCount);
	if (want_fasta)
		out_fasta.open(opt(fasta), ThreadCount);
	if (want_nuhexfasta)
		out_nuhex.open(opt(nuhexfasta), ThreadCount);
	if (want_kappafasta)
		out_kappa.open(opt(kappafasta), ThreadCount);

	s_fast_need_nu_on_chain =
		(want_can || want_nuhexfasta || want_kappafasta);
	s_fast_need_nu =
		(s_fast_need_nu_on_chain || want_bcb);

	s_fast_src = 0;
	s_fast_bca = (want_bca ? &out_bca : 0);
	s_fast_bcb = (want_bcb ? &out_bcb : 0);
	s_fast_cal = (want_cal ? &out_cal : 0);
	s_fast_can = (want_can ? &out_can : 0);
	s_fast_fasta = (want_fasta ? &out_fasta : 0);
	s_fast_nuhex = (want_nuhexfasta ? &out_nuhex : 0);
	s_fast_kappa = (want_kappafasta ? &out_kappa : 0);
	s_fast_N = 0;
	s_fast_minlen = s_MinChainLength;
	s_fast_next = 0;
	s_fast_done = 0;
	s_fast_input = 0;
	s_fast_converted = 0;
	s_fast_tooshort = 0;
	s_fast_shortest = UINT_MAX;
	s_fast_cal_in = fin;
	s_fast_cal_insz = sz;
	s_fast_cal_nthreads = ThreadCount;

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts.push_back(new thread(FastCalThreadBody, ThreadIndex));
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		ts[ThreadIndex]->join();
		delete ts[ThreadIndex];
		}
	uint done = s_fast_done.load();
	Progress("%s chains converted\n", FloatToStr(done));

	if (want_cal)
		FinalizeTextShards(out_cal, "CAL");
	if (want_can)
		FinalizeTextShards(out_can, "CAN");
	if (want_fasta)
		FinalizeTextShards(out_fasta, "FASTA");
	if (want_nuhexfasta)
		FinalizeTextShards(out_nuhex, "nuhexfasta");
	if (want_kappafasta)
		FinalizeTextShards(out_kappa, "kappafasta");
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
	CloseStdioFile(fin);
	s_fast_cal_in = 0;
	s_fast_cal_insz = 0;

	ProgressLog("\n");
	ProgressLogPrefix("%u / %u converted (%s)\n",
		s_fast_converted, s_fast_input, IntToStr(s_fast_input));
	if (s_fast_tooshort > 0)
		ProgressLogPrefix("%u too short (min length %u, shortest %u)\n",
			s_fast_tooshort, s_fast_minlen,
			(s_fast_shortest == UINT_MAX ? 0 : s_fast_shortest));

	s_fast_bca = 0;
	s_fast_bcb = 0;
	s_fast_cal = 0;
	s_fast_can = 0;
	s_fast_fasta = 0;
	s_fast_nuhex = 0;
	s_fast_kappa = 0;
	s_fast_need_nu = false;
	s_fast_need_nu_on_chain = false;
	}

static bool FastPathEligible(bool want_cal, bool want_can, bool want_bca,
	bool want_bcb, bool want_fasta, bool want_nuhexfasta,
	bool want_kappafasta, bool want_pdbcaoutdir, bool have_labels)
	{
	if (want_pdbcaoutdir)
		return false;
	if (!(want_cal || want_can || want_bca || want_bcb ||
		want_fasta || want_nuhexfasta || want_kappafasta))
		return false;
	if (have_labels)
		return false;
	if (!IsRegularFile(g_Arg1))
		return false;
	string Ext;
	GetExtFromPathName(g_Arg1, Ext);
	ToLower(Ext);
	return (Ext == "bca" || Ext == "bcb" || Ext == "cal");
	}

static void RunFastBcx(bool want_cal, bool want_can, bool want_bca,
	bool want_bcb, bool want_fasta, bool want_nuhexfasta,
	bool want_kappafasta)
	{
	BCAData src;
	src.Open(g_Arg1);
	const uint N = src.GetChainCount();

	BCAData out_bca;
	BCAData out_bcb;
	text_shards_t out_cal;
	text_shards_t out_can;
	text_shards_t out_fasta;
	text_shards_t out_nuhex;
	text_shards_t out_kappa;
	const uint ThreadCount = GetRequestedThreadCount();
	if (want_bca)
		out_bca.CreateSharded(opt(bca), false, ThreadCount);
	if (want_bcb)
		out_bcb.CreateSharded(opt(bcb), true, ThreadCount);
	if (want_cal)
		out_cal.open(opt(cal), ThreadCount);
	if (want_can)
		out_can.open(opt(can), ThreadCount);
	if (want_fasta)
		out_fasta.open(opt(fasta), ThreadCount);
	if (want_nuhexfasta)
		out_nuhex.open(opt(nuhexfasta), ThreadCount);
	if (want_kappafasta)
		out_kappa.open(opt(kappafasta), ThreadCount);

	s_fast_need_nu_on_chain =
		(want_can || want_nuhexfasta || want_kappafasta);
	s_fast_need_nu =
		(s_fast_need_nu_on_chain || want_bcb);

	s_fast_src = &src;
	s_fast_bca = (want_bca ? &out_bca : 0);
	s_fast_bcb = (want_bcb ? &out_bcb : 0);
	s_fast_cal = (want_cal ? &out_cal : 0);
	s_fast_can = (want_can ? &out_can : 0);
	s_fast_fasta = (want_fasta ? &out_fasta : 0);
	s_fast_nuhex = (want_nuhexfasta ? &out_nuhex : 0);
	s_fast_kappa = (want_kappafasta ? &out_kappa : 0);
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

	if (want_cal)
		FinalizeTextShards(out_cal, "CAL");
	if (want_can)
		FinalizeTextShards(out_can, "CAN");
	if (want_fasta)
		FinalizeTextShards(out_fasta, "FASTA");
	if (want_nuhexfasta)
		FinalizeTextShards(out_nuhex, "nuhexfasta");
	if (want_kappafasta)
		FinalizeTextShards(out_kappa, "kappafasta");
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
	s_fast_cal = 0;
	s_fast_can = 0;
	s_fast_fasta = 0;
	s_fast_nuhex = 0;
	s_fast_kappa = 0;
	s_fast_need_nu = false;
	s_fast_need_nu_on_chain = false;
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

	s_MinChainLength = flat_params::get_min_chainlength();

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
		string Ext;
		GetExtFromPathName(g_Arg1, Ext);
		ToLower(Ext);
		if (Ext == "cal")
			RunFastCal(want_cal, want_can, want_bca, s_want_bcb,
				want_fasta, want_nuhexfasta, s_want_kappafasta);
		else
			RunFastBcx(want_cal, want_can, want_bca, s_want_bcb,
				want_fasta, want_nuhexfasta, s_want_kappafasta);
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

void cmd_convert()
	{
	cmd_flat_convert();
	}

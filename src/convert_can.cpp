#include "myutils.h"
#include "flat_params.h"
#include "flat_chain_reader.h"
#include "chaq.h"
#include "seqdb.h"

extern uint8_t g_nucode_to_kappacode[256];

void cmd_convert_structs_to_can()
	{
	const string &chainfn = g_Arg1;
	if (optset_output) Die("Use -can not -output");
	if (!optset_can) Die("Must specify -can OUTPUTFILE");

	const uint maxL = 4000;
	const uint M = flat_params::m_distmx_bandwidth;
	sid_t *distmx = myalloc(sid_t, maxL*M);
	uint8_t *codeseq_nu = myalloc(uint8_t, maxL);
	uint8_t *codeseq_kappa = myalloc(uint8_t, maxL);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, maxL);

	flat_chain_reader CR;
	CR.Open(g_Arg1);

	FILE *fcan = CreateStdioFile(opt(can));

	uint nchain = 0;
	for (;;)
		{
		const flat_chain_t *chain = CR.GetNext();
		if (chain == 0) break;
		++nchain;
		if (nchain%1000 == 0) Progress("%u chains converted\r", nchain);
		uint L = chain->get_length();
		asserta(L < maxL);//TODO
		chaq::fill_codeseq_nu_from_chain(
			chain, distmx, &cv, codeseq_nu, maxL);
		const char *charseq_aa20 = chain->m_aa->m_data;
		fprintf(fcan, ">%s\n", chain->m_label.c_str());
		for (uint pos = 0; pos < L; ++pos)
			{
			char aa = chain->get_aa(pos);
			float x, y, z;
			chain->get_coords(pos, x, y, z);
			uint8_t code_nu = codeseq_nu[pos];
			fprintf(fcan, "%c\t%.1f\t%.1f\t%.1f\t%02x\n", aa, x, y, z, code_nu);
			}
		}
	Progress("%u chains converted\n", nchain);

	CloseStdioFile(fcan);
	}

void cmd_convert_can_to_kappa_fasta()
	{
	if (optset_output) Die("Use -fasta not -output");
	if (!optset_fasta) Die("Must specify -fasta FILENAME");

	FILE *f = CreateStdioFile(opt(fasta));

	string stem;
	GetStemName(g_Arg1, stem);
	Progress("Reading %s...", stem.c_str());
	vector<string> lines;
	ReadLinesFromFile(g_Arg1, lines);
	Progress(" done\n");

	const uint nline = uint(lines.size());
	string label;
	string kappa_seq;
	kappa_seq.reserve(1000);
	vector<string> flds;
	uint nseq = 0;
	for (uint lineidx = 0; lineidx < nline; ++lineidx)
		{
		const string &line = lines[lineidx];
		if (line.size() == 0)
			continue;
		if (line[0] == '>')
			{
			if (lineidx > 0)
				{
				if (label.size() == 0)
					Die("Empty label at line %u", lineidx+1);

				if (kappa_seq.size() == 0)
					Die("Zero-length structure >%s at line %u",
						label.c_str(), lineidx);
				}
			SeqToFasta(f, label, kappa_seq);
			label = line.substr(1);
			kappa_seq.clear();
			++nseq;
			if (nseq%1000 == 0)
				Progress("%u sequences\r", nseq);
			continue;
			}

		// 0          1       2       3     4
		// P       37.6    13.7    32.7    64
		Split(line, flds, '\t');
		if (flds.size() != 5)
			Die("Expected five fields in line %u, got %u",
				lineidx+1, uint(flds.size()));
		if (flds[0].size() != 1)
			Die("Expected one aa letter in line %u, got '%s'",
				lineidx+1, flds[0].c_str());
		const string &hex = flds[4];
		if (hex.size() != 2)
			Die("Expected 2-digit hex in line %u, got '%s'",
				lineidx+1, hex.c_str());
		char *endptr = 0;
		long nu = strtol(hex.c_str(), &endptr, 16);
		if (endptr == hex.c_str())
			Die("Expected hex digits in line %u, got '%s'",
				lineidx+1, hex.c_str());
		if (nu < 0 || nu >= 256)
			Die("Expected hex digits in line %u, got '%s'=%ld",
				lineidx+1, hex.c_str(), nu);
		uint8_t nu_code = uint8_t(nu);
		assert(long(nu_code) == nu);

		uint8_t kappa_code = g_nucode_to_kappacode[nu_code];
		assert(kappa_code < 32);
		char kappa_char = g_LetterToCharMu[kappa_code];
		kappa_seq += kappa_char;
		}
	Progress("%u sequences\n", nseq);
	CloseStdioFile(f);
	}

void cmd_convert_structs_to_bca()
	{
	const string &chainfn = g_Arg1;
	if (optset_output) Die("Use -bca not -output");
	if (!optset_bca) Die("Must specify -bca OUTPUTFILE");

	const uint maxL = 4000;
	const uint M = flat_params::m_distmx_bandwidth;
	sid_t *distmx = myalloc(sid_t, maxL*M);
	uint8_t *codeseq_nu = myalloc(uint8_t, maxL);
	uint8_t *codeseq_kappa = myalloc(uint8_t, maxL);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, maxL);

	flat_chain_reader CR;
	CR.Open(g_Arg1);

	BCAData BCA;
	BCA.Create(opt(bca));

	uint nchain = 0;
	for (;;)
		{
		const flat_chain_t *chain = CR.GetNext();
		if (chain == 0) break;
		++nchain;
		if (nchain%1000 == 0) Progress("%u chains read\r", nchain);
		BCA.write_flat_chain(chain);
		}
	Progress("%u chains read\n", nchain);
	Progress("finalizing... ");
	BCA.Close();
	Progress("done\n");
	}

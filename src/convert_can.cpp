#include "myutils.h"
#include "flat_params.h"
#include "flat_chain_reader.h"
#include "chaq.h"

void cmd_convert_can()
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

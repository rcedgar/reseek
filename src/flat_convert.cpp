#include "myutils.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_chain_reader.h"
#include "bcadata.h"
#include "chaq.h"

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

void cmd_flat_convert()
	{
	if (optset_output)
		Die("Use -cal, -can, -bca, -bcb, -fasta, -hexfasta or -kappafasta not -output");

	const bool want_cal = optset_cal;
	const bool want_can = optset_can;
	const bool want_bca = optset_bca;
	const bool want_bcb = optset_bcb;
	const bool want_fasta = optset_fasta;
	const bool want_hexfasta = optset_hexfasta;
	const bool want_kappafasta = optset_kappafasta;

	if (!want_cal && !want_can && !want_bca && !want_bcb &&
		!want_fasta && !want_hexfasta && !want_kappafasta)
		Die("Must specify one or more output options: "
		  "-cal, -can, -bca, -bcb, -fasta, -hexfasta, -kappafasta");

	FILE *fcal = CreateStdioFile(opt(cal));
	FILE *fcan = CreateStdioFile(opt(can));
	FILE *ffasta = CreateStdioFile(opt(fasta));
	FILE *fhexfasta = CreateStdioFile(opt(hexfasta));
	FILE *fkappafasta = CreateStdioFile(opt(kappafasta));

	BCAData BCA;
	BCAData BCB;
	if (want_bca)
		BCA.Create(opt(bca), false);
	if (want_bcb)
		BCB.Create(opt(bcb), true);

	chaq_vecs2 cv;
	uint8_t *codeseq_kappa = 0;
	if (want_bcb)
		chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
	if (want_kappafasta)
		codeseq_kappa = myalloc(uint8_t, flat_params::m_maxL);

	flat_chain_reader CR;
	CR.Open(g_Arg1);

	uint nchain = 0;
	for (;;)
		{
		const flat_chain_t *chain = CR.GetNext();
		if (chain == 0)
			break;
		++nchain;
		if (nchain%1000 == 0)
			Progress("%u chains converted\r", nchain);

		const uint L = chain->get_length();
		asserta(L < flat_params::m_maxL);
		asserta(chain->has_nu());
		const uint8_t *nu = chain->get_nu_data();

		chain->to_cal(fcal);
		WriteCan(fcan, chain);
		chain->to_fasta(ffasta);
		codeseq_to_hexfasta(fhexfasta, chain->m_label, nu, L);
		WriteKappaFasta(fkappafasta, chain, codeseq_kappa);

		if (want_bca)
			BCA.write_flat_chain(chain, &cv);
		if (want_bcb)
			BCB.write_flat_chain(chain, &cv);
		}

	Progress("%u chains converted\n", nchain);

	CloseStdioFile(fcal);
	CloseStdioFile(fcan);
	CloseStdioFile(ffasta);
	CloseStdioFile(fhexfasta);
	CloseStdioFile(fkappafasta);

	if (want_bca)
		{
		Progress("finalizing BCA... ");
		BCA.Close();
		Progress("done\n");
		}
	if (want_bcb)
		{
		Progress("finalizing BCB... ");
		BCB.Close();
		Progress("done\n");
		}

	if (want_bcb)
		chaq::free_chaq_vecs2(cv);
	myfree(codeseq_kappa);
	}

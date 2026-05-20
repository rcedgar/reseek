#include "myutils.h"
#include "statsig.h"
#include "parabench.h"
#include "flat_profiles.h"

void cmd_nubench()
	{
	string chainfn, hexfafn;
	if (EndsWith(g_Arg1, ".hexfa"))
		hexfafn = g_Arg1;
	else if (EndsWith(g_Arg1, ".bca"))
		chainfn = g_Arg1;
	else
		Die("must be .hexfa or .bca");

	flat_params params;

	Paralign::set_final_nu();
	Paralign::LogMatrix();

	ParaBench PS;
	if (optset_lookup)
		PS.ReadLookup(opt(lookup));

	flat_profiles fp;
	vector<flat_chain_t *> chains;
	if (chainfn != "")
		{
		vector<string> alpha_names;
		alpha_names.push_back("aa20");
		alpha_names.push_back("aa4");
		alpha_names.push_back("pm2");
		alpha_names.push_back("sec32");
		params.set_alpha_names(alpha_names);
		asserta(optset_lookup);
		read_flat_chains(chainfn, chains);
		fp.from_chains_lookup(*PS.m_look, chains);
		fp.set_nu_codeseqs();
		PS.m_ByteSeqs.clear();
		uint ndom = PS.m_look->get_ndom();
		PS.m_ByteSeqs.resize(ndom);
		for (uint domidx = 0; domidx < ndom; ++domidx)
			{
			const uint8_t *codeseq = fp.m_nu_codeseqs[domidx];
			const uint L = fp.get_length(domidx);
			PS.m_ByteSeqs[domidx].reserve(L);
			for (uint pos = 0; pos < L; ++pos)
				PS.m_ByteSeqs[domidx].push_back(codeseq[pos]);
			}
		}
	else if (hexfafn != "")
		PS.GetByteSeqs(hexfafn, "nuletters"); // hexfa
	else
		asserta(false);

	if (!optset_lookup)
		PS.SetLookupFromLabels();
	PS.Search("para", false);
	PS.SetScoreOrder();
	PS.Bench();
	PS.WriteHits(opt(output), opt(include_self), opt(triangle));
	}

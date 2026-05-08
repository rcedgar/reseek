#include "myutils.h"
#include "chaq.h"
#include "flat_chain.h"
#include "alpha.h"
#include "quantize.h"

void cmd_flat_quantized2fa()
	{
	Die("TODO");
	asserta(optset_fasta);
	asserta(optset_alpha_size);
	asserta(optset_feature);
	const uint alpha_size = opt(alpha_size);
	const string &chainfn = g_Arg1;
	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);
	FILE *f = CreateStdioFile(opt(fasta));
	const string &feature = opt(feature);
	FAN fan = str2FAN(feature.c_str());

	asserta(is_quantized(fan));
	const uint16_t *thresholds = chaq::get_thresholds(fan, alpha_size);
	const uint16_t undef_value = chaq::get_undef_value(fan, alpha_size);

	for (uint i = 0; i < nchain; ++i)
		{
		const flat_chain_t *chain = chains[i];
		uint L = chain->get_length();
		char *charseq = myalloc(char, L);
		chaq::slow_get_charseq_binned(chain, fan, alpha_size,
				thresholds, undef_value, charseq);
		SeqToFasta(f, chains[i]->m_label.c_str(), charseq, L);
		myfree(charseq);
		}
	CloseStdioFile(f);
	}

void cmd_flat_discrete2fa()
	{
	asserta(optset_fasta);
	asserta(optset_alpha_size);
	asserta(optset_feature);
	const uint alpha_size = opt(alpha_size);
	const string &chainfn = g_Arg1;
	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);
	FILE *f = CreateStdioFile(opt(fasta));
	const string &feature = opt(feature);
	FAN fan = str2FAN(feature.c_str());

	asserta(!is_quantized(fan));
	uint8_t undef_code = chaq::get_undef_code(fan, alpha_size);

	for (uint i = 0; i < nchain; ++i)
		{
		const flat_chain_t *chain = chains[i];
		uint L = chain->get_length();
		char *charseq = myalloc(char, L);
		chaq::slow_get_charseq_discrete(chain, fan, alpha_size,
				undef_code, charseq);
		SeqToFasta(f, chains[i]->m_label.c_str(), charseq, L);
		myfree(charseq);
		}
	CloseStdioFile(f);
	}

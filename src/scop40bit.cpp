#include "myutils.h"
#include "scop40bench.h"
#include "pdbchain.h"
#include "sort.h"

/***
==> 3dblastaln <==
d12asa_ d12asa_ 1981
d12asa_ d1nnha_ 428
d12asa_ d1eova2 406

==> 3dblastswaln <==
                    <<< blank line!!  
d12asa_ d12asa_ 1981
d12asa_ d1nnha_ 428

==> cealn <==
d12asa_ d12asa_ 8.03 0.00 1.00
d12asa_ d1nnha_ 6.81 2.40 0.81
d12asa_ d1b8aa2 6.70 2.76 0.75

==> cleswaln <==
d12asa_ d12asa_ 15284
d12asa_ d1b8aa2 5739
d12asa_ d1nnha_ 5620

==> dalialn <==
d12asa_ d12asa_ 58.0 0.0
d12asa_ d1nnha_ 27.0 2.7
d12asa_ d1b8aa2 25.1 2.7

==> foldseekaln <==
      0       1     2     3   4 5   6     7 8     9        10     11
d1a1xa_	d1a1xa_	0.000	106	105	0	1	106	1	106	2.549E-163	1045
d1a1xa_	d1jsga_	0.000	108	105	0	1	106	4	111	9.137E-93	622
d1a1xa_	d3saoa_	0.000	77	76	0	26	102	27	103	2.448E-21	188

==> mmseqsaln <==
d12asa_ d12asa_ 674
d12asa_ d2gz4a1 39
d12asa_ d1b8aa2 30

==> tmaln <==
d12asa_ d12asa_ 1.0
d12asa_ d1b8aa2 0.75576
d12asa_ d1nnha_ 0.751037

==> tmfastaln <==
d12asa_ d12asa_ 1.0000 1.0000
d12asa_ d1b8aa2 0.7558 0.7394
d12asa_ d1nnha_ 0.7510 0.8309
***/

void cmd_scop40bit()
	{
	const string &Algo = g_Arg1;
	asserta(optset_input);
	asserta(optset_output);

	SCOP40Bench SB;
	asserta(optset_benchmode);
	SB.m_Mode = string(opt_benchmode);
	SB.ReadChains(opt_input);
	SB.BuildDomFamIndexesFromQueryChainLabels();
	SB.ReadHits_Tsv(Algo);
	SB.WriteBit(opt_output);
	}

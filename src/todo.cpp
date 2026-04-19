/***
@@TODOs

reverse profiles not chains (handedness)

DALI score
	Optimize by gapless hsp pairs in nu d.p. mx, greedily build?
	Transition to DALI score for low mega, better for fold recognition
	Weight DALI score by SS/SEC alphabet match

TM-score, how to do without Kabsch?

FEN etc

PUU

Use landmarks to make a composition vector
	|----------|------------|----------|
	 ^^^^^^^^^^
	Frequency vector
	of Kappa or Nu
	=> dot product "BLOSUM" score for fast
	domain identification, align to top few

RNA structures

Turn off exceptions
	GCC/Clang: -fno-exceptions (and often -fno-rtti if desired)
	MSVC: /EHs-c- (disable C++ exceptions)

Cache Mu in db binary for prefilter kmers

Test scripts and data for very short chains e.g. <10aa.

Use pair-wise alignment to anchor contact map profile alignment.

Test speedup with __restrict for flat_base::m_data

Position-specific gap penalties.

Redefine neighbors after local alignment is constructed.

Low-complexity weighting

>d1g9ga_/a.102.1.2 (nendist)
PJAAPJAAPPAAPPAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAA

Feature Plot Structure Aligner (FPSA)
https://docs.google.com/document/d/1YRq1LQcEIgraHt_5G_V5-YU1YIwEmnu_MHpQQZmQsuY/edit?tab=t.0#heading=h.pcvo5gmuhifv
***/
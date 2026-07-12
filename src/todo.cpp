/***
TODOs

Remove ref-counting.
Exploit cached nu.

reverse profiles not chains (handedness)

Classify pair as fam / sf / fold, re-align
	this is good if it improves accuracy by all truth standards.

DALI score
	Optimize by gapless hsp pairs in nu d.p. mx, greedily build?
	Transition to DALI score for low mega, better for fold recognition
	Weight DALI score by SS/SEC alphabet match

Fold recognition by image recognition of Reseek similarity matrix (NO)

Fold recognition by "tableaux":
	Consensus fold represented as "contact map" of SSEs (inc. loop)
	  separated by linkers
	Each segment of the structure has feature vector
		Frequency of sec32 letters
		mean & stdev length
	Each cell in the map has relation e.g.
		non-local contact (for pair of SSEs)
		if adjacent then { linear, hairpin, kink... }
	
	Learn consensus from MStA

TM-score, how to do without Kabsch?
	rotfreetm.cpp/h
	binning longer distances possibly improves speed of refinement
	  (but not calculation given alignment)
	re-define alignment without 1:1, e.g. residue-> 3D 
	  coordinate of chain defined by C-alpha links.

FEN etc

PUU domain detection

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

turnd and angle have almost same information?

Low-complexity weighting

>d1g9ga_/a.102.1.2 (nendist)
PJAAPJAAPPAAPPAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAAPAAA

Feature Plot Structure Aligner (FPSA)
https://docs.google.com/document/d/1YRq1LQcEIgraHt_5G_V5-YU1YIwEmnu_MHpQQZmQsuY/edit?tab=t.0#heading=h.pcvo5gmuhifv
***/
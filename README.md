<p align="left"><img src="https://drive5.com/images/reseek_logo.jpg" height="100"/></p>

Reseek is a novel protein structure alignment algorithm which improves sensitivity in protein homolog detection
compared to state-of-the-art methods including DALI, TM-align and Foldseek with improved speed over Foldseek, the
fastest previous method.

Reseek is based on sequence alignment where each residue in the protein backbone is represented by a 
letter in a novel “mega-alphabet” of 85,899,345,920 (∼10<sup>11</sup>) distinct states.

Method sensitivity was measured on the SCOP40 benchmark using superfamily as the truth standard, focusing
on the regime with false-positive error rates <10 per query, corresponding to E<10 for an ideal E-value.

This is a preview beta release, new features and improved documentation will hopefully follow soon.
Feedback is welcome via github Issues.

<pre>
All-vs-all alignment (excluding self-hits)
    reseek -search STRUCTS -mode MODE -output hits.tsv 

Search query structures against database
    reseek -search Q_STRUCTS -db DB_STRUCTS -mode MODE -output hits.tsv

Align two structures
    reseek -search NAME1.pdb -db NAME2.pdb -mode MODE -aln aln.txt

Output options for -search
   -aln FILE     # Alignments in human-readable format
   -output FILE  # Hits in tabbed text format with 8 fields:
                 #   Evalue Query Target
                 # (More output formats coming soon)

Search and alignment options
  -mode MODE     # veryfast|fast|sensitive (default fast)
  -evalue E      # Max E-value (default report all alignments)
  -omega X       # Omega accelerator (floating-point)
  -minu U        # K-mer accelerator (integer)
  -gapopen X     # Gap-open penalty (floating-point >= 0, default 1.1)
  -gapext X      # Gap-extend penalty (floating-point >= 0, default 0.14)
  -dbsize D      # Effective database size for E-value (default actual size)
  -usort         # U-sort accelerator (default off)
  -maxaccepts N  # If U-sort, max hits <= E-value (default 1)
  -maxrejects N  # If U-sort, max hits > E-value (default 32)

Convert PDB file(s) to .cal (C-alpha) format
    reseek -pdb2cal STRUCTS -output structs.cal

STRUCTS argument is one of:
   NAME.pdb      # PDB file (mmCIF support will be added soon)
   NAME.files    # Text file with PDB file/pathnames, one per line
   NAME.cal      # C-alpha (.cal) file, recommended for databases
</pre>


### Reference

Edgar, Robert C. "Sequence alignment using large protein structure alphabets improves sensitivity to remote homologs" [https://www.biorxiv.org/content/10.1101/2024.05.24.595840v1](https://www.biorxiv.org/content/10.1101/2024.05.24.595840v1)


### SCOP40 benchmark code and results

https://github.com/rcedgar/reseek_scop40

<p align="center"><img src="https://github.com/rcedgar/reseek_scop40/blob/master/results/sens_vs_err_sf.png" align="left" width="700"/></p>
<p></p>

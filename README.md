![Reseek](http://drive5.com/images/reseek_logo.jpg)

Reseek is a novel protein structure alignment algorithm which doubles sensitivity in protein homolog detection
compared to state-of-the-art methods including DALI, TM-align and Foldseek with improved speed over Foldseek, the
fastest previous method. 
Reseek is based on sequence alignment where each residue in the protein backbone is represented by a 
letter in a novel “mega-alphabet” of 85,899,345,920 (∼10<sup>11</sup>) distinct states.

This is a preview beta release, new features and improved documentation will hopefully follow soon.
Feedback is welcome via github Issues.

<pre>
All-vs-all alignment (excluding self-hits)
    reseek -search STRUCTS -mode MODE -output hits.txt 

Search query against database
    reseek -search Q_STRUCTS -db DB_STRUCTS -mode MODE -output hits.txt

Align two structures
    reseek -search NAME1.pdb -db NAME2.pdb -mode MODE -aln aln.txt

Output options for -search
   -aln FILE     # Alignments in human-readable format
   -output FILE  # Hits in tabbed text format with 8 fields:
                 #   Evalue Query Target Qstart Qend Tstart Tend CIGAR
                 # (More output formats coming soon)

Search and alignment options
  -mode MODE     # veryfast|fast|sensitive|verysensitive (required)
  -evalue E      # Max E-value (default 10)
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
Edgar, Robert C. "Sequence alignment using large protein structure alphabets doubles sensitivity to remote homologs"     
[https://www.biorxiv.org/content/10.1101/2024.05.24.595840v1](https://www.biorxiv.org/content/10.1101/2024.05.24.595840v1)

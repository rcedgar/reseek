![Reseek](http://drive5.com/images/reseek_logo2.jpg)

Reseek is a protein structure search and alignment algorithm which improves sensitivity in protein homolog detection
compared to state-of-the-art methods including [DALI](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3749), [TM-align](https://academic.oup.com/nar/article-abstract/33/7/2302/2401364) and [Foldseek](https://www.biorxiv.org/content/10.1101/2022.02.07.479398.abstract). Speed is often similar to Foldseek while memory requirements are often much less, though scaling trade-offs are different.

### Online structure search

Search a protein structure against AFDB, PDB or BFVD with typical results in 2 to 5 minutes.

<hr>

[https://reseek.online](https://reseek.online)

<hr>

### What's new in version 3

Search accuracy is substantially improved over previous releases, and search is faster with large databases.

The underlying statistical model is enhanced, now offering three sets of parameters optimized for SCOP family, superfamily and fold relationships, respectively

**Families** are somewhat arbitrary groups with similar functions and clear amino acid sequence similarity.

**Superfamilies** are groups of families which are probably homologous.

**Folds** are groups of families with similar tertiary structure.

For more information about these categories, see [this SCOP documentation page](https://www.ebi.ac.uk/pdbe/scop/about).

The `-stats MODEL` option is required in version 3, `MODEL` is `family`, `superfamily` or `fold`. For example, the family model has a larger weight for amino acid sequence similarity, while the fold model has higher weights for tertiary structure features.

The preferred measure of statistical significance is P-value, specified by the `pvalue` column. E-values are deprecated because they are misleading for structure searches. P-values are adjusted according to a null model tuned to the chosen truth standard, i.e. family, superfamily or fold.

### Reseek achieves highest accuracy in homolog detection and statistical significance estimate

On the [SCOP40 benchmark test](https://www.pnas.org/doi/abs/10.1073/pnas.95.11.6073) (see results later below), Reseek has substantially higher ability to discriminate homologs compared to previous algorithms including [DALI](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3749), [TM-align](https://academic.oup.com/nar/article-abstract/33/7/2302/2401364) and [Foldseek](https://www.biorxiv.org/content/10.1101/2022.02.07.479398.abstract).

Reseek also provides a more accurate estimate of statistical significance, enabling users to set a cutoff based on an acceptable number of false positives for a given search, while Foldseek E-values may over-estimate significance by 5 to 6 orders of magnitude (reference below).

Paper is here: [https://drive5.com/reseek/Reseek3-2026-07-16_preprint.pdf](https://drive5.com/reseek/Reseek3-2026-07-16_preprint.pdf).

### YouTube talk describing the algorithm

Reseek is based on sequence alignment where each residue in the protein backbone is represented by a letter in a novel “mega-alphabet” of 85,899,345,920 (∼10<sup>11</sup>) distinct structure states. This talk explains how it works.

[<img src="https://drive5.com/reseek/youtube_snip.gif" width="150">](https://www.youtube.com/watch?v=BzIgqdm9xDs)

### Command line
<pre>
Commands
  -search        # Alignment (e.g. DB search, pairwise, all-vs-all)
  -convert       # Convert file formats (e.g. create DB)
  -alignpair     # Pair-wise alignment and superposition

Search against database
    reseek -search STRUCTS -db db.bcb -output hits.tsv
                 # use -convert to create .bcb files (below)

Align and superpose two structures
    reseek -alignpair 1XYZ.pdb -input2 2ABC.pdb
           -aln FILE     # Sequence alignment (text)
           -output FILE  # Rotated 1XYZ (PDB format)

Output options for -search
   -aln FILE     # Alignments in human-readable format
   -output FILE  # Hits in tabbed text format
   -columns name1+name2+name3...
                 # Output columns, names are:
                 #   query   Query label
                 #   target  Target label
                 #   qlo     Start of aligment in query
                 #   qhi     End of aligment in query
                 #   tlo     Start of aligment in target
                 #   thi     End of aligment in target
                 #   ql      Query length
                 #   tl      Target length
                 #   pctid   Percent identity of alignment
                 #   cigar   CIGAR string
                 #   pvalue  P-value
                 #   qrow    Aligned query sequence with gaps (local)
                 #   trow    Aligned target sequence with gaps (local)
                 #   qrowg   Aligned query sequence with gaps (global)
                 #   trowg   Aligned target sequence with gaps (global)
                 #   qcovpct Query coverage (percent)
                 #   tcovpct Target coverage (percent)
                 #   std     query+target+qlo+qhi+ql+tlo+thi+tl+pctid+pvalue (default)

Search and alignment options
  -fast or -sensitive  # Required
  -pvalue P            # Max P-value (default 0.001)

Convert between file formats
    reseek -convert STRUCTS [one or more output options]
           -cal FILENAME    # .cal format, text with a.a. and C-alpha x,y,z
           -bca FILENAME    # .bca format, binary .cal
           -bcb FILENAME    # .bca format, binary .cal with nu, required for db
           -fasta FILENAME  # FASTA format

STRUCTS argument is one of:
   NAME.cif or NAME.mmcif     # PDBx/mmCIF file
   NAME.pdb                   # Legacy format PDB file
   NAME.cal                   # C-alpha tabbed text format with chain(s)
   NAME.bcb                   # Binary C-alpha with nu, required for db
   NAME.bca                   # Binary C-alpha (used by earlier reseek versions)
   NAME.files                 # Text file with one STRUCT per line,
                              #   may be filename, directory or .files
   DIRECTORYNAME              # Directory (and its sub-directories) is searched
                              #   for known file types including .pdb, .files etc.

Other options:
   -log FILENAME              # Log file with errors, warnings, time and memory.
   -threads N                 # Number of threads, default number of CPU cores.

More documentation at https://drive5.com/reseek
</pre>

#### Build from source on Linux x86
<pre>
cd src/; chmod +x build_linux_x86.bash ; ./build_linux_x86.bash
</pre>

#### Build from source on Windows
Load `reseek.vcxproj` into Microsoft Visual Studio and use the Build command.

#### Ignore static link warning
Don't worry about a warning something like this, it's expected:
<pre>
warning: Using 'dlopen' in statically linked applications requires
  at runtime the shared libraries from the glibc version used for linking
</pre>

### More documentation

[https://drive5.com/reseek](https://drive5.com/reseek)

### Searching very large databases
Reseek and Foldseek achieve orders of magnitude faster searches of large database such as AFDB compared to previous state-of-the-art methods DALI and TM-align. In terms of accuracy, my benchmark results show that Reseek-fast (i.e., with the `-fast`) option has higher accuracy than all other algorithms by most metrics, but scaling compared to Foldseek involves quite different trade-offs (preprint in preparation). Foldseek indexes the whole database in memory, which sometimes enables sub-linear scaling with number of query structures but requires a large amount of RAM. Reseek indexes the query, thereby using much less RAM with linear scaling in query size. With more than ~1,000 query structures, for Reseek it is recommended to split the query into subsets of ~1k and search these separately. The results below illustrate the point. Random subsets of the PDB with 100, 1,100 and 10,000 structures were searched against a subset of AFDB which was clustered at 50% aa similarity. Experiments were run on a machine with 780Gb RAM using 64 threads. Results are shown below. With 100 structures, reseek is faster, but with 1,000+ structures Foldseek is faster.

<pre>
  query   ____________Elapsed time_______________
structs   foldseek  reseek-fast  reseek-sensitive
    100         9m          6m                23m
   1000        19m         61m             4h:45m  
  10000      2h33m         12h                60h

  query   _____________Max memory________________
structs   foldseek  reseek-fast  reseek-sensitive
    100      227Gb        4.4Gb             4.5Gb
   1000      284Gb        6.6Gb             6.7Gb
  10000      289Gb            -                 -
</pre>

### Database downloads

[abdb50.bcb](https://serratus-public.s3.us-east-1.amazonaws.com/rce/reseek_dbs/afdb50.bcb) &nbsp;&nbsp;AFBD50 AFDB clustered at 50% aa identity 53.7M structures (107Gb)  

[esmdb30.bcb](https://serratus-public.s3.us-east-1.amazonaws.com/rce/reseek_dbs/esmdb30.bcb) &nbsp;&nbsp;ESMDB30 ESM Altas clustered at 30% aa identity 34.6M structures (55Gb)  

[bfvd.bcb](https://serratus-public.s3.us-east-1.amazonaws.com/rce/reseek_dbs/bfvd.bcb) &nbsp;&nbsp;Big Fine Virus DB 347k structures (486Mb)  

[pdb.bcb](https://serratus-public.s3.us-east-1.amazonaws.com/rce/reseek_dbs/pdb.bcb) &nbsp;&nbsp;PDB 900k structures (1.7Gb)

[scop40.bcb](https://serratus-public.s3.us-east-1.amazonaws.com/rce/reseek_dbs/scop40.bcb) &nbsp;&nbsp;SCOP40 11,211 structures (16Mb)

[scop40x.bcb](https://serratus-public.s3.us-east-1.amazonaws.com/rce/reseek_dbs/scop40x.bcb) &nbsp;&nbsp;Curated SCOP40 8,291 structures (13Mb)

[cath40.bcb](https://serratus-public.s3.us-east-1.amazonaws.com/rce/reseek_dbs/cath40.bcb) &nbsp;&nbsp;CATH40 34,647 structures (41Mb)

### SCOP40 benchmark code and results
https://github.com/rcedgar/reseek_bench

![Reseek](https://drive5.com/images/reseek3_accuracy_preprint_fig.jpg)

**Accuracy plots for SCOP40c superfamily and fold.**   
The figure gives CVE, PR and ROC plots for the tested algorithm
using SCOP40c (curated SCOP40, [https://github.com/rcedgar/scop40c](https://github.com/rcedgar/scop40c)) as reference and superfamily and fold as truth standards. For a CVE plot, a lower curve is better while
for PR and ROC a higher curve is better. These curves show that Reseek-sensitive has higher accuracy than other
tested methods in the high-scoring regime (up to 10 errors per query for CVE, ≥ 70% precision for PR and < 10−4
FPR for ROC).

### References

Edgar RC. "Protein structure alignment by Reseek improves sensitivity to remote homologs" (_Bioinformatics_ 2024) Nov;40(11):btae687. 
[https://academic.oup.com/bioinformatics/article/40/11/btae687/7901215](https://academic.oup.com/bioinformatics/article/40/11/btae687/7901215)

Edgar RC. and Sahakyan S. "Protein structure alignment significance is often exaggerated" (_bioRxiv_ 2025) [https://www.biorxiv.org/content/10.1101/2025.07.17.665375v1](https://www.biorxiv.org/content/10.1101/2025.07.17.665375v1)

Edgar RC. Rich structure alphabets enable highest accuracy protein search (_bioRxiv_ 2026) [https://www.biorxiv.org/cgi/content/short/2026.07.24.740611v1](https://www.biorxiv.org/cgi/content/short/2026.07.24.740611v1)

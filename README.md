# bamrefine

[![PyPI version](https://badge.fury.io/py/bamrefine.svg)](https://badge.fury.io/py/bamrefine)

This is a BAM filtering/masking program that  masks the positions from an
alignment if that read carries a variant position that could be affected
by post-mortem damage (PMD). e.g. If a read carries a C/\* position at its
5' end or it carries a G/\* position at its 3' end. The regions treated as 5'
or 3' "end" is not pre-defined and should be provided to the program with its
argument (see `-l` or `--pmd-length-threshold` parameter)

It is designed to be an alternative to hard clipping first N bases 
of BAM reads with the suspicion of those regions containing 
non-authentic (post mortem damage artifact) SNPs, or 
ignoring transition SNPs alltogether. Thus, avoiding a substantial amount
of data loss that is encountered when using mainstream pipelines.

## Install

The easiest way to `bamrefine` is from PyPI using `pip`:

```
pip install bamrefine
```

you can also install `bamrefine` from a local download:

to do so, first download the latest release from the [github releases page](https://github.com/etkayapar/bamRefine/tags), and uncompress the downloaded archive. Once you have the source code downloaded, you can navigate into the folder and use pip to install from the copy you have downloaded:

```
cd bamRefine-0.1.X
pip install .
```


Upon installation of the package, the command `bamrefine` should be available to you.

## Masking strategy

Consider the following example.

![](https://github.com/etkayapar/bamRefine/blob/master/figs/healthyAlignment_scaled.jpg?raw=true)

above, in the left, you see the table of 
a theoretical SNP catalog with the 
variant positions and minor and major alleles
observed in a population. In the right,you see a theoretical alignment 
of three paired end reads to the reference (orange) region covering all
the variants listed on the left. Some reads match the reference, others
carry variants (green).

Below, a more realistic example for ancient DNA reads, in
addition to the expected biological variation among the aligned
reads, there is also unexpedted variants resulting from DNA damage
at the ends of fragmented ancient DNA. 

![](https://github.com/etkayapar/bamRefine/blob/master/figs/damagedAlignment_scaled.jpg?raw=true)

At the 5' ends, these damages will occur as C-\>T changes, and at the 
3' ends they will appear as G-\>A changes. 

![](https://github.com/etkayapar/bamRefine/blob/master/figs/PMD_smiley.jpg?raw=true)


Therefore, a C/\* variant in the table becomes a possible PMD source
if the read overlaps with it at the reads 5' end, and this is also
the case for a G/\* variant if a read overlaps with it at the 3' end.

Bamrefine, first characterizes the input SNP collection into two tables
as 5' end suspects and 3' end suspects:

![](https://github.com/etkayapar/bamRefine/blob/master/figs/snpTables.jpg?raw=true) 

and then iterates over the reads in the BAM file searching for overlaps
with the two suspect tables at the appropriate ends of the reads (overlaps
with the 5' end of the read and the 5' table variants, 3' end of the reads 
and the 3' table variants). Then if a read carries a suspected/dangerous variant 
at the ends, the program masks those specific position in the read so that
they won't interfere with the SNP calling for the downstream steps:

![](https://github.com/etkayapar/bamRefine/blob/master/figs/maskingExample.jpg?raw=true)


As you can see at the right hand side example of the above picture, 
there is still one unexpected 'A' (in red), it stayed unmasked even 
though it is within the 3' end lookup range (indicated with
double-headed arrow). This is because, 

- there is no entry in the
  SNP collection at that position of the genome. Meaning, since there
  won't be any SNP calling concerning that specific locus, there is no
  need to bother with unexpected variants like this.
- `bamrefine` does not check the specific allele that a read
  carries, rather does the masking based on the alignment overlap of the
  reads ends and the "dangerous" positions in the SNP catalog. See the
  5' end maskings at the  left part of above example: Even though the read
  carries 'C' bases (not different from the reference), those 'C's are
  still masked.

This therefore means:

- reads that do not overlap with the problematic variant 
  positions at their ends will still provide information 
  regarding the allele pool
- you disregard unreliable information from
  the alignments while preserving anything
  else. i.e. If there is more than one variant spanning
  a single read and one is a possible PMD source, by 
  selectively masking the problematic position, you still
  have the reliable information coming from the other 
  position.

## Usage

```bamrefine [parameters] [flags] <in.bam> <out.bam>```

### Parameters:

  * `-s, --snps`: SNP collection file. This can be either 4/5 column BED 
  (genomic position [with or without the `start` position column] + 
  Minor and Major allele) or SNP (that is used by EIGENSTRAT) formatted files. 
  Format distinction is done by checking the file extension so it is important 
  that the input file follows the format of the extension it has. See the 
  `sample_data` data directory in the installation path for examples of each of 
  the three supported file types. You can find the installation directory by running
  `pip3 show bamrefine` or `pip show bamrefine`.
  * `-p, --threads`: Threads to run the program in parallel.
  * `-l, --pmd-length-threshold`: Either a single integer (e.g `-l 10`) or two integers 
  separated by a comma (e.g `-l 2,0`) representing different length values corresponding
  5' and 3' ends , respectively. Positions that are up to and including this far in any read will 
  be evaluated and masked.
  Also, it is recommended that you have PMD related statistics (e.g. up to which position there is a high 
  risk of seeing a PMD artifact) regarding your libraries before you run this program on the BAM files 
  because the `--pmd-length-threshold` parameter requires an user specified value for the program to run, 
  there is no default value. This can be achieved by using [PMDtools](https://github.com/pontussk/PMDtools)


### Flags:

  * `-v, --verbose`: verbose output of progress.
  * `-t, --add-tags`: Add tags to reads in output BAM file related to masking statistics 
    using optional SAM fields in alignment records. e.g. `ZC:Z:2,1`  and `ZP:Z:0,5;-3` 
    would mean that the program masked n=2 5' and n=1 3' positions and they were at index 
    0,5 and -3 in the read sequence. 3' masking positions are represented with negative 
	indices that start counting from the 3' end of the read and is compatible with python 
	list indices.
  * `-k, --keep-tmp`: Don't remove the temprorary run directory. It will include 
    intermediate BAM files, cached SNPs, etc. this directory is under the same
    directory with specified ouptut BAM, named as `.YYYY-MM-DD_HH-MM-SS_<out.bam>_tmp_bamrefine`

  

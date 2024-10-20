<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [bamrefine](#bamrefine)
    - [Installation](#installation)
        - [Via PyPI](#via-pypi)
        - [From source](#from-source)
    - [Masking strategy](#masking-strategy)
    - [Usage](#usage)
        - [Parameters:](#parameters)
        - [Flags:](#flags)

<!-- markdown-toc end -->
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

## Installation

Regardless of the installation method you choose below, we recommend installing `bamrefine` inside a Python virtual environment as opposed to installing it at the system or user level. Depending on your environment you may actually have to use this approach if you don't have the right permissions.

We recommend you use the python module `venv` for this purpose, it is a part of The Python Standard Library thus should not require you to install it separately. You can read more about how the module works from [here](https://docs.python.org/3/library/venv.html).

First create and activate your `venv` virtual environment:

```
python3 -m venv my_env
source my_env/bin/activate
```

after this step you may proceed with either one of the below installation methods.

### Via PyPI

The easiest way to install `bamrefine` is from PyPI using `pip`. 


You can run the below line in your terminal to install `bamrefine` inside the virtual environment from PyPI using `pip`:

```
pip install bamrefine
```

Upon installation of the package, the command `bamrefine` should be available to you.

### From source

You can also install `bamrefine` from a local download of this repository.

To do so, first download the latest release (or the specific version you want to install) from the [github releases page](https://github.com/etkayapar/bamRefine/tags), and uncompress the downloaded archive.

Then, install `bamrefine` inside the virtual environment from the source code using `pip`:

```
cd bamRefine-0.2.X
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

  * `-S, --single-stranded`: Run the program in single-stranded mode. This flag should be 
    turned on for BAM files that are generated from single-stranded libraries. When turned 
    on, C/* variants are considered for overlapping read positions at both 5' and 3' ends.
    (Instead of using C/* for 5' ends and G/* for 3' ends)
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



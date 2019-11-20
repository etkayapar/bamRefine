# bamrefine

## What is it?

This is a BAM filtering/masking program that conditionally masks
transition SNPs from BAM reads. 

It designed to be an alternative to hard clipping first N bases 
of BAM reads with the suspicion of those regions containing 
non-authentic (post mortem damage artifact) SNPs, or 
ignoring transition SNPs alltogether. Instead, bamrefine masks 
only the transition SNPs in the possible PMD regions from both 
ends of BAM reads (scope of the region specified by the user)

## Usage

```bamrefine [parameters] [flags] <in.bam> <out.bam>```

### Parameters:

  * `-s, --snps`: SNP collection file. This can be either
    4/5 column BED or SNP formatted files. Fields/columns
    can be either tab separated or space delimited. Both
    are acceptable. For safety, dont mix and match two
    delimiters in the same file. For examples, see
    ```data/```
  * `-g, --genome-fasta`: Path to genome fasta file. This
    is needed only for fetcing contig names. Genome should be
    indexed beforehand with `samtools faidx`.
  * `-l, --pmd-length-threshold`: N nucleotide region from
    both ends of a read to be treated as a possible PMD region.
  * `-p, --threads`: Threads to run the program in parallel.

### Flags:

  * `-v, --verbose`: verbose output of progress.
  
## Dependencies

This is a small program written in python3, so it requires pyton3.5
or above. It has several python libraries that it depends
on:

  * `cython`
  * `pysam`

Though this is not a dependency by any means, it is recommended that you
have PMD related statistics (e.g. up to which position there is a high 
risk of seeing a PMD artifact) regarding your libraries before you run this
program on the BAM files because the `--pmd-length-threshold` parameter
requires an user specified value for the program to run. This can be achieved
by using [PMDtools](https://github.com/pontussk/PMDtools)

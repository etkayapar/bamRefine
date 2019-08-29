# bamrefine

## What is it?

This is a BAM filtering/masking program that conditionally masks
transition SNPs from BAM reads. 

It designed to be an alternative to hard clipping first N bases 
of BAM reads with the suspicion of those regions containing 
non-authentic (post mortem damage artifact) SNPs, or 
ignoring transition SNPs alltogether. Instead, bamrefine masks 
the transition SNPs in the possible PMD regions from both ends 
of BAM reads (scope of the region specified by the user)

## Usage

```bamrefine [parameters] [flags] <in.bam> <out.bam>```

Parameters:

  * ```-s, --snps```: SNP collection file. This can be either
    BED or SNP formatted files. To see an example file, see
    ```data/```
  * ```-g, --genome-fasta```: Path to genome fasta file. This
    is needed only for fetcing contig names. Genome should be
    indexed beforehand.
  * ```-l, --pmd-length-threshold```: N nucleotide region from
    both ends of a read to be treated as a possible PMD region.
  * ```-p, --threads```: Threads to run the program in parallel.

Flags:

  * ```-v, --verbose```: verbose output of progress.

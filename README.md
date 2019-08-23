# bamrefine

## What is it?

This is a BAM filtering/masking program that conditionally masks
transition SNPs from BAM reads. In a nutshell, If a SNP is a
transition and if its located on the ends of a BAM read. This
ensures that, you are not using non-authentic SNPs that are
most likely result of PMD(Post-mortem Damage). Therefore it is
recommended that you use this program having PMD scores of your
library and tweak the ```--pmd-length-threshold``` parameter 
accordingly.

## Usage

bamrefine [parameters] [flags] <in.bam> <out.bam>

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

  * ``-v, --verbose```: verbose output of progress.

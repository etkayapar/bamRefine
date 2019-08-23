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

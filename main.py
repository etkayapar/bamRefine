from bamRefine import *
import pysam

## Pseudo ------
for each snp in snps:
    get pileups for that pos, if there is a variant and it is
    at the ends of the read, mask the read, or else skip the read.


def snpMask(snp, samF):

    crm = snp[0]
    loc = int(snp[1]) - 1
    ref = snp[2]
    alt = snp[3]

    for p_column in samF.pileup(crm, loc, loc+1):
        for p_read in p_column.pileups:
            p_alignment = p_read.alignment
            if p_alignment.query_sequence[p_read.query_position] == alt and
            p_read.query_position:
                continue

















import sys
from bamRefine_cy import *
# import pysam


inName = sys.argv[1]
contig = sys.argv[2]



snps = parseSNPs('chr_pos_ref_alt_1240K.all.snp')
# snps = list(snps.values())
# snps = snps[:3000]

ouName = contig + ".bam"

inBAM = pysam.AlignmentFile(inName, 'rb')
ouBAM = pysam.AlignmentFile(ouName, 'wb', template=inBAM)

processBAM(inBAM, ouBAM, snps, contig)

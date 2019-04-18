import sys
from bamRefine import *
import pysam


inName = sys.argv[1]
contig = sys.argv[2]



snps = parseSNPs('chr_pos_ref_alt_1240K.all.snp')
# snps = list(snps.values())
# snps = snps[:3000]

ouName = contig + ".bam"

inBAM = pysam.AlignmentFile(inName, 'rb')
ouBAM = pysam.AlignmentFile(ouName, 'wb', template=inBAM)

#fReads = snpMask(snps, inBAM)

#filterBAM(inBAM, ouBAM, fReads)


for read in inBAM.fetch(contig):
    bamL = read.to_dict()
    mask, m_pos = flagReads(snps, bamL)
    if mask == True:
        for p in m_pos:
            # print('in read %s, position %d' % (bamL['name'], p))
            t = list(bamL['seq']) ; t[p] = 'N' ; t = "".join(t)
            bamL['seq'] = t
        # bamLine --> pysam.AlignedSegment
        bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
        ouBAM.write(bamL)
    else:
        bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
        ouBAM.write(bamL)

inBAM.close()
ouBAM.close()

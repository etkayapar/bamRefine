from bamRefine import *
import pysam
import dill
from pathos.multiprocessing import ProcessingPool

snps = parseSNPs('chr_pos_ref_alt_1240K.all.snp')
# snps = list(snps.values())
# snps = snps[:3000]

inBAM = pysam.AlignmentFile('../deneme.bam', 'rb')
ouBAM = pysam.AlignmentFile('../denemeout.bam', 'wb', template=inBAM)

#fReads = snpMask(snps, inBAM)

#filterBAM(inBAM, ouBAM, fReads)


#for read in inBAM.fetch():




def filterBAMline(read):
    global snps
    global inBAM
    # bamL = read.to_dict()
    bamL = read
    mask, m_pos = flagReads(snps, bamL)
    if mask == True:
        for p in m_pos:
            #print('in read %s, position %d' % (bamL['name'], p))
            t = list(bamL['seq']) ; t[p] = 'N' ; t = "".join(t)
            bamL['seq'] = t
        # bamLine --> pysam.AlignedSegment
        #bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
        #ouBAM.write(bamL)
        return bamL
    else:
        #bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
        #ouBAM.write(bamL)
        return bamL

# bamF = [read.to_dict() for read in inBAM.fetch()]

with ProcessingPool(4) as p:
    results = p.map(filterBAMline, inBAM.fetch() , chunksize = 50)

print('Finished')


# with open('maskedReads.txt', 'w') as f:
#     for r in fReads.keys():
#         f.write(r + "\n")


print('writing bam')
for read in results:
    bamL = pysam.AlignedSegment.from_dict(read, inBAM.header)
    ouBAM.write(bamL)


# res = list(results)


inBAM.close()
ouBAM.close()

from bamRefine.functions import *
#import os

dirN = os.path.dirname(os.path.realpath(__file__))
os.chdir(dirN)



snps = "random_reads.snp"
snps = handleSNPs(snps)
lookup_l = 10
lookup_r = 10

def test_flagReads_complex():
    inName = "complex_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps, bamL, lookup_l, lookup_r, read)

        assert mask == "mask" 
        assert m_pos == [1,5,8,-9]
        assert m_side == [0,0,0,1]
    inBAM.close()

def test_flagReads_complex_asymmetric():
    inName = "complex_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps, bamL, 10, 8, read)

        assert mask == "mask" 
        assert m_pos == [1,5,8]
        assert m_side == [0,0,0]
    inBAM.close()

def test_flagReads_short():
    inName = "short_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps, bamL, lookup_l, lookup_r, read)
        print(m_pos)
        print(m_side)
        print(snps)

        assert mask == "mask" 
        assert m_pos == [0,1,5,-2,-5]
        assert m_side == [0,0,0,1,1]
    inBAM.close()


    

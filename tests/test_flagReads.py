from bamRefine.functions import *
import os

dirN = os.path.dirname(os.path.realpath(__file__))
os.chdir(dirN)



lookup_l = 10
lookup_r = 10

### Double-stranded mode (default) -------------------------
snps_ds = "random_reads.snp"
snps_ds = parseSNPs(snps_ds, False)
def test_flagReads_complex():
    inName = "complex_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps_ds, bamL, lookup_l, lookup_r, read)

        assert mask == "mask" 
        assert m_pos == [1,5,8,-9]
        assert m_side == [0,0,0,1]
    inBAM.close()

def test_flagReads_complex_asymmetric():
    inName = "complex_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps_ds, bamL, 10, 8, read)

        assert mask == "mask" 
        assert m_pos == [1,5,8]
        assert m_side == [0,0,0]
    inBAM.close()

def test_flagReads_short():
    inName = "short_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps_ds, bamL, lookup_l, lookup_r, read)
        print(m_pos)
        print(m_side)
        print(snps_ds)

        assert mask == "mask" 
        assert m_pos == [0,1,5,-2,-5]
        assert m_side == [0,0,0,1,1]
    inBAM.close()
## ---------------------------------------------------

## single-stranded mode
snps_ss = "random_reads_ss.snp"
snps_ss = parseSNPs(snps_ss, singleStranded=True)

def test_flagReads_complex_single():
    inName = "complex_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps_ss, bamL, lookup_l, lookup_r, read)

        assert mask == "mask" 
        assert m_pos == [1,5,8]
        assert m_side == [0,0,0]
    inBAM.close()

def test_flagReads_short_single():
    inName = "short_read.sam"
    inBAM = pysam.AlignmentFile(inName, 'rb')
    for read in inBAM.fetch():
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps_ss, bamL, lookup_l, lookup_r, read)
        print(m_pos)
        print(m_side)
        print(snps_ds)

        assert mask == "mask" 
        assert m_pos == [0,1,5,-5,-9,-10]
        assert m_side == [0,0,0,1,1,1]
    inBAM.close()

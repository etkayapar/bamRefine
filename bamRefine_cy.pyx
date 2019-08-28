import pysam
import os
import pickle
import sys
import subprocess as sp

cdef isTransition(str ref,  str bamR):

    table = {
        'A': 'G',
        'T': 'C',
        'G': 'A',
        'C': 'T'}

    try:
        if bamR == table[ref]:
            return True
        else:
            return False
    except KeyError:
        return False

def parseBam(bamL, fields = ('chrm', 'pos', 'seq')):
    bamL = bamL.strip().split()
    allF = {
        'qname': bamL[0].encode('utf-8'),
        'flag':  int(bamL[1]),
        'chrm':  bamL[2].encode('utf-8'),
        'pos':   int(bamL[3]),
        'mapq':  int(bamL[4]),
        'cigar': bamL[5].encode('utf-8'),
        'seq':   bamL[9].encode('utf-8'),
        'qual':  bamL[10].encode('utf-8')}

    # return tuple([allF[f] for f in fields])
    return allF

# Parse snp catalogue

cpdef flagReads(snpLocDic, bamLine, look):

    '''
    Find if there is a snp in the loaded SNP catalogue
    for the given range. It is usally for checking if a
    SAM/BAM read contains a SNP.
    '''

    # chrm, start, seq = bamLine
    cdef str chrm = bamLine['ref_name']
    cdef int start = int(bamLine['ref_pos'])
    cdef str seq = bamLine['seq']
    cdef int end = start + len(seq) #end variable is 1 more then actual end position for ease of use in range() calls
    cdef int ref = 2
    cdef int alt = 3
    cdef int nt
    #cdef str key
    cdef list snp
    cdef list inspectRange

    cdef list snpList = [] # store transitions pos. in the ends

    if len(seq) > look:
        inspectRange = list(range(start, start+look)) + list(range(end-look-1, end))
    else:
        inspectRange = list(range(start, end))

    #sys.displayhook(inspectRange)
    for nt in inspectRange:
        key = chrm + " " + str(nt)
        try:
            snp = snpLocDic[key]
            if snp[ref] == seq[nt-start]:
                continue
            else:
                snpList.append(nt-start)

        except KeyError:
            continue

    if len(snpList) > 0:
        return ('mask', snpList)
    else:
        return ('nomask', snpList)

def parseSNPs(fName):
    snpF = open(fName)
    # curC = 'chr1'
    # chrm = 'chr1'
    snps = {}
    if fName.endswith('.snp'):
        chrI = 1
        posI = 3
        refI, altI = (4, 5)
    elif fName.endswith('.bed'):
        ncol_cmdList = ["cat", fName, "|", "wc -w"]
        ncol_cmd = " ".join(ncol_cmdList)
        p = sp.Popen([ncol_cmd], shell=True,
                     stdout = sp.PIPE,
                     stderr = sp.STDOUT)
        ncol, _ = p.communicate()
        ncol = int(ncol.decode('utf-8').strip())

        chrI, posI, refI, altI = tuple(range(4))

        if ncol == 5:
           posI, refI, altI = tuple([x+1 for x in (posI, refI, altI)]) 


    for snp in snpF:
        snp = snp.strip().split()
        curC = snp[chrI]
        if not isTransition(snp[refI], snp[altI]):
            # ignoring transversions
            continue
        key = curC + " " + snp[posI]
        snps[key] = [snp[x] for x in [chrI, posI, refI, altI]]
    snpF.close()
    return snps


def parseSNPs_old(fName):
    snpF = open(fName)
    # curC = 'chr1'
    # chrm = 'chr1'
    chrm = '1'
    tmpS = []
    snps = {}

    for snp in snpF:
        snp = snp.strip().split()
        curC = snp[0]
        if not isTransition(snp[2], snp[3]):
            # ignoring transversions
            continue
        try:
            snps[curC].append(snp[1:])
        except KeyError:
            snps[curC] = []
            snps[curC].append(snp[1:])

    snpF.close()
    return snps



def snpMask(snps, samF):

    faultyReads = {}
    for snp in snps:
        crm = snp[0]
        loc = int(snp[1]) - 1
        ref = snp[2]
        alt = snp[3]

        for p_column in samF.pileup(crm, loc, loc+1):
            if p_column.nsegments > 0:
                for p_read in p_column.pileups:
                    if not p_read.is_refskip and not p_read.is_del:
                        p_alignment = p_read.alignment
                        readLen = len(p_alignment.query_sequence)
                        lookup = list(range(10)) + list(range(readLen-11, readLen))
                        pos = p_read.query_position
                        var = p_alignment.query_sequence[p_read.query_position] == alt
                        if var and (pos in lookup):
                            aln = p_alignment
                            masked = aln.query_sequence.replace(str(pos),'N')
                            aln.query_sequence = masked
                            faultyReads[p_alignment.query_name] = aln

    return faultyReads

def filterBAM(inAln, outAln, faultyReads):

    for read in inAln.fetch():
        try:
            outAln.write(faultyReads[read.query_name])
        except KeyError:
            outAln.write(read)


def processBAM(inBAM, ouBAM, snps, contig, lookup):

    # cdef dict bamL
    cdef list m_pos
    cdef list t
    cdef str t1
    cdef list q
    cdef str q1
    lookup = int(lookup)

    for read in inBAM.fetch(contig):
        bamL = read.to_dict()
        #print(bamL)
        mask, m_pos = flagReads(snps, bamL, lookup)
        if mask == 'mask':
            for p in m_pos:
                # print('in read %s, position %d' % (bamL['name'], p))
                t = list(bamL['seq']) ; t[p] = 'N' ; t1 = "".join(t) ; del(t)
                bamL['seq'] = t1
                q = list(bamL['qual']) ; q[p] = '!' ; q1 = "".join(q) ; del(q)
                bamL['qual'] = q1
                # bamLine --> pysam.AlignedSegment
            bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
            ouBAM.write(bamL)
        elif mask == 'nomask':
            bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
            ouBAM.write(bamL)
        else:
            continue

    inBAM.close()
    ouBAM.close()


def handleSNPs(fName):
    pickleName = os.path.basename(fName) + ".brf"
    f_exists = os.path.isfile(pickleName)
    if f_exists:
        f = open(pickleName, 'rb')
        snps = pickle.load(f)
        f.close()
    else:
        snps = parseSNPs(fName)
        f = open(pickleName, 'wb')
        pickle.dump(snps, f)
        f.close()

    return snps

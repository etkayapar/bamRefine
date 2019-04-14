def isTransition(ref, bamR):

    # table = {
    #     'A': ['G', 'X'],
    #     'T': ['C', 'X'],
    #     'G': ['A', 'X'],
    #     'C': ['T', 'X'],
    #     'X': ['A', 'T', 'G', 'C', 'X']}

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
        'qname': bamL[0],
        'flag':  int(bamL[1]),
        'chrm':  bamL[2],
        'pos':   int(bamL[3]),
        'mapq':  int(bamL[4]),
        'cigar': bamL[5],
        'seq':   bamL[9],
        'qual':  bamL[10]}

    # return tuple([allF[f] for f in fields])
    return allF

# Parse snp catalogue

def flagReads(snpLocDic, bamLine):

    '''
    Find if there is a snp in the loaded SNP catalogue
    for the given range. It is usally for checking if a 
    SAM/BAM read contains a SNP.
    '''

    # chrm, start, seq = bamLine
    chrm = bamLine['ref_name']
    start = int(bamLine['ref_pos'])
    seq = bamLine['seq']
    end = start + len(seq)
    ref = 2
    alt = 3
    lookup = list(range(start, start+10)) + list(range(end-9, end))

    snpList = [] # store transitions pos. in the ends

    for nt in lookup:
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
        return (True, snpList)
    else:
        return (False, snpList)

def parseSNPs(fName):
    snpF = open(fName)
    # curC = 'chr1'
    # chrm = 'chr1'
    snps = {}

    for snp in snpF:
        snp = snp.strip().split()
        curC = snp[0]
        if not isTransition(snp[2], snp[3]):
            # ignoring transversions
            continue
        key = curC + " " + snp[1]
        snps[key] = snp
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














from sys import stdin


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
        'qname' = bamL[0],
        'flag' = int(bamL[1]),
        'chrm' = bamL[2],
        'pos' = int(bamL[3]),
        'mapq'= int(bamline[4]),
        'cigar' = bamL[5],
        'seq' = bamL[9],
        'qual' = bamL[10]}

    return tuple([allF[f] for f in fields])

# Parse snp catalogue

def flagReads(snpLocDic, bamLine):

    '''
    Find if there is a snp in the loaded SNP catalogue
    for the given range. It is usally for checking if a 
    SAM/BAM read contains a SNP.
    '''

    chrm, start, seq = bamLine
    end = start + len(seq)
    ref = 1
    alt = 2
    lookup = list(range(start, start+11)) + list(range(end-10, end+1))

    snpList = [] # store transitions pos. in the ends

    for nt in lookup:
        key = chrm + " " + str(nt)
        try:
            snp = snpLocDic[key]
            if snp[ref] == seq[nt-start-1]:
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
        snps[key] = snp[1:]
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


snps = parseSNPs('chr_pos_ref_alt_1240K.all.snp')

with open('../deneme.sam') as sam:
    for i,l in enumerate(sam):
        msg = flagReads(snps, parseBam(l))
        if msg[0]:
            p = str(i) + '\t' + '\t'.join([str(x) for x in msg[1]])
            print(p)



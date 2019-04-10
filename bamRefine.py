from sys import stdin


def isTransition(ref, bamR):

    table = {
        'a': 'g',
        't': 'c',
        'g': 'a',
        'c': 't'}

    if table[ref] == bamR:
        return True
    else:
        return False

def parseBam(bamLine):
    bamLine = bamLine.strip().split()
    pos = bamLine[3]
    seq = bamLine[9]

    return (bamLine, pos, seq)

# Parse snp catalogue

def parseSNPs(fName):
    snpF = open(fName)
    # curC = 'chr1'
    chrm = 'chr1'
    tmpS = []
    snps = {}

    for snp in snpF:
        snp = snp.strip().split()
        curC = snp[0]
        # print("current: %s, chrm: %s" % (curC, chrm))
        try:
            snps[curC].append(snp[1:])
        except KeyError:
            snps[curC] = []
            snps[curC].append(snp[1:])

    snpF.close()
    return snps





# for bamL in stdin:
#     bamL, seq, pos = parseBam(bamL)
#     # Find if seq includes SNPs


from sys import stdin

# Parse snp catalogue

def parseSNP_locs(fName):
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
            snps[curC].append(int(snp[1]))
        except KeyError:
            snps[curC] = []
            snps[curC].append(int(snp[1]))

    snpF.close()
    return snps

def dists(snpLocDic):
    chrms = snpLocDic.keys()
    retval = [] 
    for c in chrms:
        tmpDists = []
        for i in range(1, len(snpLocDic[c])):
            dist = snpLocDic[c][i] - snpLocDic[c][i-1]
            tmpDists.append(dist)
        retval += tmpDists

    return retval


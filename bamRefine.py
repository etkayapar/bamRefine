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



for bamL in stdin:
    bamL, seq, pos = parseBam(bamL)
    # Find if seq includes SNPs


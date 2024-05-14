import pysam
import os
import pickle
import sys
import subprocess as sp
from importlib.metadata import version

def isDangerous(var):
    if ("C" in var) and ("G" in var):
        return (True, 'both')
    elif ("C" in var):
        return (True, '0')
    elif ("G" in var):
        return (True, '1')
    else:
        return (False, None)

def isDangerous_single(var):
    if ("C" in var):
        return (True, 'both')
    else:
        return (False, None)

def generateTags(snpList, sideList):
    tag1 = ",".join([str(sideList.count(x)) for x in [0,1]])
    tag1 = ('ZC',tag1)

    sideSwitchI = "".join(str(x) for x in sideList).find('1')
    if sideSwitchI != -1:
        tag2_1 = ",".join([str(x) for x in snpList[0:sideSwitchI]])
        tag2_2 = ",".join([str(x) for x in snpList[sideSwitchI:]])
    else:
        tag2_1 = ",".join([str(x) for x in snpList])
        tag2_2 = ""
    tag2 = tag2_1 + ";" + tag2_2
    tag2 = ('ZP', tag2)

    return [tag1, tag2]

def modifyHeader(input_header, bamrefine_version, bamrefine_command_line):
    output_header_dict = input_header.to_dict()
    prev_bamrefine = [x['ID'] for x in output_header_dict['PG'] if 'bamrefine' in x['ID']]
    this_bamrefine = ""
    if len(prev_bamrefine) > 0:
        if '.' in prev_bamrefine[-1]:
            last_bamrefine = prev_bamrefine[-1]
            this_bamrefine = last_bamrefine.split(".")[-1]
            this_bamrefine = str(int(this_bamrefine) + 1)
            this_bamrefine += "."
        else:
            this_bamrefine = ".1"
    ID = "bamrefine" + this_bamrefine

    output_header_dict['PG'].append({'ID': ID,
                                     'PN': 'bamrefine',
                                     'VN': bamrefine_version,
                                     'PP': output_header_dict['PG'][-1]['ID'],
                                     'CL': bamrefine_command_line})

    output_header = pysam.AlignmentHeader.from_dict(output_header_dict)

    return output_header

def flagReads(snpLocDic, bamLine, look_l, look_r, bamRecord):

    '''
    Flag positions to be masked. Mapped read positions in <bamRecord> that
    overlap with selected variants in the <snpLocDic> within <look_l> bases
    from the 5' end and <look_r> bases from the 3' end will be flagged for
    masking.
    '''

    chrm = bamLine['ref_name']
    seq = bamLine['seq']
    snpList = [] # store positions to be masked
    sideList= [] # store mask sides of positions (5' or 3')

    ## make look_l = look_r if the latter not set:
    look_r = {True: look_l, False: look_r}[look_r is None]
    look_list = [look_l, look_r]

    refpos = bamRecord.get_reference_positions(full_length=True)
    refpos = [x + 1 if x is not None else x for x in refpos] ## pysam uses 0-based indices
    read_len = len(seq)
    ## Lookup range should be about the physical first N bases in the read, not the first
    ## N reference bases. I will assume this until I finish implementing this feature.
    ## this might change later.

    ## Default inspectRange:
    inspectRange = [refpos[:look_l], refpos[read_len-1:read_len-look_r-1:-1]]

    ## Adjust if read is too short for the lookup values from either side:
    for side in range(2):
        sign = [1,-1][side]
        if read_len <= look_list[side]:
            inspectRange[side] = refpos[::sign] ## reverse the pos list if 3' end

    for side in range(2):
        for i, nt in enumerate(inspectRange[side]):
            shift = [0,1][side]
            sign = [1,-1][side]
            i = i + shift
            i = i * sign ### get indices for the 3' end
            key = chrm + " " + str(nt)
            try:
                snp = snpLocDic[side * -1][key] ## Always use the 5' side if ss-mode
            except KeyError:
                continue
            else:
                snpList.append(i)
                sideList.append(side)

    if len(snpList) > 0:
        return ('mask', snpList, sideList)
    else:
        return ('nomask', snpList, sideList)

def parseSNPs(fName):
    snpF = open(fName)
    # curC = 'chr1'
    # chrm = 'chr1'
    snps = [{}, {}]
    if fName.endswith('.snp'):
        chrI = 1
        posI = 3
        refI, altI = (4, 5)
    elif fName.endswith('.bed'):
        ncol_cmdList = ["head -n1", fName, "|", "wc -w"]
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
        dangerous, side = isDangerous([snp[refI], snp[altI]])
        if not dangerous:
            # ignoring otherwise
            continue
        key = curC + " " + snp[posI]
        if side != 'both':
            snps[int(side)][key] = [snp[x] for x in [chrI, posI, refI, altI]]
        else:
            snps[0][key] = [snp[x] for x in [chrI, posI, refI, altI]]
            snps[1][key] = [snp[x] for x in [chrI, posI, refI, altI]]

    snpF.close()
    return snps

def processBAM(inBAM, ouBAM, snps, contig, lookup, addTags = False):

    lookup = lookup.split(",")
    lookup_l = int(lookup[0])
    try:
        lookup_r = int(lookup[1])
    except IndexError:
        lookup_r = None

    statsF = open(contig+"_stats.txt", 'w')
    stats = [0,0]

    for read in inBAM.fetch(contig):
        bamL = read.to_dict()
        mask, m_pos, m_side = flagReads(snps, bamL, lookup_l, lookup_r, read)
        if mask == 'mask':
            t = list(bamL['seq'])
            q = list(bamL['qual'])
            for p in m_pos:
                # print('in read %s, position %d' % (bamL['name'], p))
                t[p] = 'N'
                q[p] = '!'
            t1 = "".join(t)
            q1 = "".join(q)
            bamL['seq'] = t1
            bamL['qual'] = q1
            bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
            st = [m_side.count(x) for x in [0,1]]
            for x in range(2):
                stats[x] += st[x]
            if addTags:
                bamL.tags += generateTags(m_pos, m_side)
            ouBAM.write(bamL)
        elif mask == 'nomask':
            bamL = pysam.AlignedSegment.from_dict(bamL, inBAM.header)
            ouBAM.write(bamL)
        else:
            continue

    ## statsF.write("5p_total," + "3p_total"+ '\n')
    statsF.write(str(stats[0]) +  "\t" + str(stats[1]) + '\n')
    inBAM.close()
    ouBAM.close()
    statsF.close()


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


def createBypassBED(inName, chrms, snps):
    snps = handleSNPs(snps)
    s_chrms = set()

    for i in range(2):
        for snp in snps[i]:
            chrm = snp.split()[0]
            if chrm in chrms:
                s_chrms.add(snp.split()[0])


    toBypass = set(chrms).difference(s_chrms)
    toFilter = s_chrms
    if toBypass != set():
        bed = pysam.idxstats(inName).split("\n")
        bed = [x.split("\t") for x in bed][:-2]
        bed = [x[0:2] for x in bed if x[2] != '0']
        bed = [x for x in bed if x[0] in toBypass]
        bed = ["\t".join([x[0], "0", x[1]]) for x in bed]
        with open("bypass.bed", 'w') as bedF:
            for line in bed:
                bedF.write(line + '\n')
        bypass = True
    else:
        toFilter = chrms
        bypass = False

    return (list(toFilter), bypass)

def fetchChromosomes(inName):
    idstats = pysam.idxstats(inName).split("\n")
    idstats = [x.split("\t") for x in idstats][:-2] ## last two contigs are meaningless
    chrms = [x[0] for x in idstats if x[2] != '0']
    return chrms

def main(args=None):

    inName = sys.argv[1]
    contig = sys.argv[2]
    lookup = sys.argv[3]
    snps = sys.argv[4]
    addTags = bool(int(sys.argv[5]))
    singleStranded = bool(int(sys.argv[6]))


    snps = handleSNPs(snps)

    ouName = contig + ".bam"

    bamrefine_commandline = os.environ['BAMREFINE_CMDLINE']
    bamrefine_version = version("bamrefine")

    inBAM = pysam.AlignmentFile(inName, 'rb')
    output_header = modifyHeader(inBAM.header, bamrefine_version, bamrefine_commandline)

    ouBAM = pysam.AlignmentFile(ouName, 'wb', header=output_header)

    processBAM(inBAM, ouBAM, snps, contig, lookup, addTags)

    return 0

"""
author: Carrie Wessinger

this generates input files for lepmap3 linkage map software

assumes vcf file is filtered for HQ biallelic snps

"""

def generatePedigree(nsamples, nF2s, crossname, grandpa, grandma, pa, ma, F2names):
    '''writes the pedigree lines'''
    # line 1 lists family name for all samples
    postfile.write('CHR\tPOS')
    for i in range(nsamples + 2):
        postfile.write('\t' + crossname)
    postfile.write('\n')
    # line 2 lists name of all samples
    postfile.write('CHR\tPOS\t' + grandpa + '\t' + grandma + '\t' + pa + '\t' + ma)
    for i in F2names:
        postfile.write('\t' + i)
    postfile.write('\n')
    # line 3 lists "father" of each individual. (is 0 for grandparents)
    postfile.write('CHR\tPOS\t0\t0\t' + grandpa + '\t' + grandpa)
    for i in range(nF2s):
        postfile.write('\t' + pa)
    postfile.write('\n')
    # line 4 lists "mother" of each individual. (is 0 for grandparents)
    postfile.write('CHR\tPOS\t0\t0\t' + grandma + '\t' + grandma)
    for i in range(nF2s):
        postfile.write('\t' + ma)
    postfile.write('\n')
    # line 5 lists "sex" of each individual. (is 0 means unknown)
    postfile.write('CHR\tPOS\t1\t2\t1\t2')
    for i in range(nF2s):
        postfile.write('\t0')
    postfile.write('\n')
    # line 6 lists "phenotype" of each individual. (not currently used)
    postfile.write('CHR\tPOS')
    for i in range(nsamples + 2):
        postfile.write('\t0')
    postfile.write('\n')
    return


def skipHeader(line):
    """ For a line in VCF file, if line is not a header line, extract columns into list """
    cols = line.replace('\n','').split('\t')
    if len(cols) < 2:
        return 'header'
    elif cols[0] == '#CHROM':
        return 'header'
    else:
        return cols

def extractVCFfields(sampleData):
    """ Extract data from sample-level fields in VCF file """
    fields = sampleData.split(':')
    if fields[0] != './.':
        phreds = fields[4].split(',')
        return phreds
    else:
        return 'missing'

def genoLs(phreds, refAllele, altAllele):
    """ assign genotype likelihoods
        Order is LLs = ['AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT']
    """
    AA = '0'
    AC = '0'
    AG = '0'
    AT = '0'
    CC = '0'
    CG = '0'
    CT = '0'
    GG = '0'
    GT = '0'
    TT = '0'
    if phreds[0] == '0':
        LRR = '1'
    else:
        LRR = '{:0.3e}'.format(10 ** (-float(phreds[0])/10.0))
    if phreds[1] == '0':
        LRA = '1'
    else:
        LRA = '{:0.3e}'.format(10 ** (-float(phreds[1])/10.0))
    if phreds[2] == '0':
        LAA = '1'
    else:
        LAA = '{:0.3e}'.format(10 ** (-float(phreds[2])/10.0))
    if refAllele == 'A':
        AA = LRR
        if altAllele  == 'C':
            AC = LRA
            CC = LAA
        elif altAllele == 'G':
            AG = LRA
            GG = LAA
        elif altAllele == 'T':
            AT = LRA
            TT = LAA
    elif refAllele == 'C':
        CC = LRR
        if altAllele == 'A':
            AC = LRA
            AA = LAA
        elif altAllele == 'G':
            CG = LRA
            GG = LAA
        elif altAllele == 'T':
            CT = LRA
            TT = LAA
    elif refAllele == 'G':
        GG = LRR
        if altAllele == 'A':
            AG = LRA
            AA = LAA
        elif altAllele == 'C':
            CG = LRA
            CC = LAA
        elif altAllele == 'T':
            GT = LRA
            TT = LAA
    elif refAllele == 'T':
        TT = LRR
        if altAllele == 'A':
            AT = LRA
            AA = LAA
        elif altAllele == 'C':
            CT = LRA
            CC = LAA
        elif altAllele == 'G':
            GT = LRA
            GG = LAA
    LLs = AA + ' ' + AC + ' ' + AG + ' ' + AT + ' ' + CC + ' ' + CG + ' ' + CT + ' ' + GG + ' ' + GT + ' ' + TT
    return(LLs)

vcf = open('allPheno_100individ.best.hwe.vcf', 'rU') # your HWE filtered vcf
F2namefile = open('all_phenoed_F2s.txt', 'rU') # a text file listing the names of just the F2s
postfile = open('allPheno_100individ_Q30.best.lepmap.input.txt', 'w') # your output file


nF2s = 386
nsamples = nF2s + 2
crossName = 'NB'
grandpa = 'bar'
grandma = 'neo'
pa = 'F1m'
ma = 'F1f'
grandpa_location = nsamples - 2 # the 0-indexed position of bar (need to add 9 to find place in vcf)
grandma_location = nsamples - 1 # the 0-indexed position of neo (need to add 9 to find place in vcf)

# READ IN NAMES OF F2s
F2names = []
for line in F2namefile:
    F2names.append(line.replace('\n', ''))

# GENERATE PEDIGREE FILE
generatePedigree(nsamples, nF2s, crossName, grandpa, grandma, pa, ma, F2names)

# GENERATE GENOTYPE FILE
for line in vcf:
    cols = line.replace('\n', '').split('\t')
    contig = cols[0]
    position = cols[1]
    ref_allele = cols[3]
    alt_allele = cols[4]
    postfile.write(contig + '\t' + position)
    colcounter = 0

    # read in data for the two grandparents
    for j in [grandpa_location, grandma_location]:
        phreds = extractVCFfields(cols[j + 9])
        if phreds != 'missing':
            LLs = genoLs(phreds, ref_allele, alt_allele)
        else:
            LLs = '1 1 1 1 1 1 1 1 1 1'
        postfile.write('\t' + LLs)
        colcounter += 1

    # create missing data for two F1 parents
    for j in range(2):
        LLs = '1 1 1 1 1 1 1 1 1 1'
        postfile.write('\t' + LLs)
        colcounter += 1

    # read in data for each F2 individual
    for j in range(nF2s):
        phreds = extractVCFfields(cols[j + 9])
        if phreds != 'missing':
            LLs = genoLs(phreds, ref_allele, alt_allele)
        else:
            LLs = '1 1 1 1 1 1 1 1 1 1'
        postfile.write('\t' + LLs)
        colcounter += 1
    postfile.write('\n')
postfile.close()

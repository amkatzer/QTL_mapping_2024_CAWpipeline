'''
author: Carrie Wessinger

this finds site coverage per sample

input is vcf file and list of samples
'''

import numpy as np

##################################
def skipHeader(line):
    cols = line.replace('\n','').split('\t')
    if len(cols) < 2:
        return 'header'
    elif cols[0] == '#CHROM':
        return 'header'
    else:
        return cols

def extractVCFfields(sampleData):
    """ Extract data from sample-level fields in VCF file """
    if sampleData.split(':')[0] != './.':
        fields = sampleData.split(':')
        if (len(fields) >= 4):
            alleleDepths = fields[1].split(',')
            totDepth = fields[2]
            phreds = fields[4].split(',')
            return [totDepth, alleleDepths, phreds]
        else:
            return 'missing'
    else:
        return 'missing'

##################################

vcf = open('20240523_CAWfilt_1individ_NB.vcf', 'rU')
outfile = open('20240523_Q30_1individ.sample.coverage.txt', 'w')
labelfile = open('all_phenoed.txt', 'rU')

nindivs = 388

###
# read in labels for samples
# store in list = "labs"

labs = []
for line in labelfile:
    cols = line.replace('\n', '').split('\t')
    labs.append([cols[0], 0, [], [], 0])  # ['name', nsites, [genos_per_site], [depth_per_site], nhetsites]

lines = 0
for line in vcf:
    cols = skipHeader(line)

    if cols != 'header':
        lines += 1
        ngenos = 0
        for j in range(9, nindivs+9):
            if cols[j] != './.':
                ngenos += 1

        for j in range(9, nindivs+9):
            annot = extractVCFfields(cols[j])
            if annot != 'missing':
                depth = annot[0]
                labs[j-9][1] += 1
                labs[j-9][2].append(ngenos)
                labs[j-9][3].append(float(depth))
                if annot[2][1] == '0':
                    labs[j-9][4] += 1
                else:
                    pass
            else:
                pass

outfile.write('sample\tsites\tmedian_genos_persite\tmedian_depth_persite\thet_freq\n')
for l in labs:
    outfile.write(l[0] + '\t' + str(l[1]) + '\t' + str(np.median(l[2])) + '\t' + str(np.median(l[3])) + '\t' + str(float(l[4])/float(l[1])) + '\n')
outfile.close()

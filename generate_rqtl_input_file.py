"""
author: Carrie Wessinger

this generates input files for rqtl linkage map software

from output files of lep-map3

"""

from numpy import mean

samplefile = open('all_phenoed_F2s.txt', 'rU')
outfile = open('Allchromos.20240524_Q30.manually.ordered.genotypes.csv', 'w')

F2s = 386
LGs = 8

markernames = []
LGnumbers = []
mappositions = []
genotypes = [[] for i in range(F2s)]

samplenames = []

for line in samplefile:
    cols = line.replace('\n', '').split('\t')
    if cols[0] != 'sampleName':
        samplenames.append(cols[0])


for i in range(LGs):
    thisfile = open('chr{0}.manually.ordered.with.genos.mapped.txt'.format(i+1), 'rU')

    for line in thisfile:
        cols = line.replace('\n', '').split('\t')

        contig = cols[0]
        position = cols[1]
        markernames.append(contig + '_' + position)

        LG = cols[2]
        LGnumbers.append(LG)

        map_pos = round(mean([float(cols[3]), float(cols[4])]), 4)
        mappositions.append(str(map_pos))

        for j in range(F2s):
            data = cols[j+5]

            if data == '1 1':
                outgeno = 'A'
            elif data == '2 2':
                outgeno = 'B'
            elif data == '1 2' or data == '2 1':
                outgeno = 'H'
            else:
                data == '?'
                print('wtf')

            genotypes[j].append(outgeno)

    thisfile.close()

# first line of output is marker names
outfile.write('id')
for m in markernames:
    outfile.write(',' + m)
outfile.write('\n')

# second line of output is linkage group numbers
for l in LGnumbers:
    outfile.write(',' + l)
outfile.write('\n')

# third line of output is map positions
for p in mappositions:
    outfile.write(',' + p)
outfile.write('\n')

# then following lines are genotype data
for i in range(F2s):
    outfile.write(samplenames[i])
    for g in genotypes[i]:
        outfile.write(',' + g)
    outfile.write('\n')

outfile.close()

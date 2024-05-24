in2a = open('allPheno_100individ_Q30.best.ids.txt', 'rU')
vcf  =open("20240523_CAWfilt_NB.vcf", "rU")
outfile = open("allPheno_100individ_Q30.best.vcf", 'w')

bestlist = []
for line in in2a:
    cols = line.replace('\n', '').split('\t')
    bestlist.append([cols[0], cols[1]])

for line in vcf:
    cols = line.replace('\n', '').split('\t')
    if len(cols) > 2 and cols[0][0] != '#':
        tig = cols[0]
        pos = cols[1]
        for x in bestlist:
            if tig == x[0] and pos == x[1]:
                outfile.write(line)
    else:
        pass
outfile.close()


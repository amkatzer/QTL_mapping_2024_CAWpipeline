# This program thins to best snp per radtag
# written by JKK

# modified by Carrie
# assumes you've already filtered for quality
# modified for use with bcftools

MinMinor = 8  # min individuals that will have minor allele

vcf = open("20240523_CAWfilt_NB.vcf", "rU")
out2a = open('allPheno_100individ_Q30.best.ids.txt', 'w')

plants = 386 + 2

last_scaff = ''
bestsnp = ''
lastpos = 0
cp = 0

bestlist = []

for line_idx, line in enumerate(vcf):
    cols = line.replace('\n', '').split('\t')
    if len(cols) > 2 and cols[0] != "#CHROM":
        # scaffold_10	6853	.	G	A	567.83	.	AC=6;AF=0.025;AN=240;BaseQRankSum=-0.678;DP=554;Dels=0.00;FS=0.000;HaplotypeScore=0.2216;InbreedingCoeff=0.0014;MLEAC=6;MLEAF=0.025;MQ=60.00;MQ0=0;MQRankSum=-1.740;QD=16.22;ReadPosRankSum=0.712	GT:AD:DP:GQ:PL	0/0:12,0:12:36:0,36,499	./.	./.	./.	./.	0/1:2,9:12:53:297,0,53	./.	0/0:6,0:6:18:0,18,246	./.	0/0:4,0:4:12:0,12,165	./.	./.	./.	0/0:4,0:4:12:0,12,166	0/0:3,0:3:9:0,9,124	0/0:7,0:7:21:0,21,283	./.	0/0:4,0:4:12:0,12,165	./.	./.	0/0:1,0:1:3:0,3,42	./.	0/0:1,0:1:3:0,3,41	./.	./.	1/1:0,3:3:9:114,9,0	./.
        scaff = cols[0]
        position = int(cols[1])
        ref_base = cols[3]
        alt_base = cols[4]
        if line_idx % 10000 == 0:
            print(scaff, position)

        mincc = 0

        if len(alt_base) == 1:  # print "why are there multiple bases here?"
            datums = [0, 0, 0]
            for j in range(9, 9 + plants):
                if cols[j].split(":")[0] != "./.":
                    geno = cols[j].split(":")[0]
                    if geno == "0/0":
                        datums[0] += 1
                    elif geno == "0/1":
                        datums[1] += 1
                    elif geno == "1/1":
                        datums[2] += 1
                    else:
                        print("wtf2 ", geno)
            mincc = min(datums[0] + datums[1], datums[2] + datums[1])

            if scaff != last_scaff or (position - lastpos) > 150:
                if cp >= MinMinor:
                    out2a.write(bestsnp)
                cp = mincc
                bestsnp = scaff + '\t' + cols[1] + '\n'
                last_scaff = scaff
                lastpos = position
            elif mincc > cp and position != lastpos:
                cp = mincc
                bestsnp = scaff + '\t' + cols[1] + '\n'

if cp >= MinMinor:
    out2a.write(bestsnp)  # last snp

out2a.close()




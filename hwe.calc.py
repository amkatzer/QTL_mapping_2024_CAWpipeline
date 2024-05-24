# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 15:04:53 2017

@author: carrie

getting basic info from vcf files

Assumes a filtered vcf!
"""

from numpy import median
from numpy import mean
from scipy.stats import chi2

in1 = open('allPheno_100individ_Q30.best.vcf', 'rU')
stats = open('allPheno_100individ_Q30.best.info.hwe.txt', 'w')

nF2s = 386
nsamples = nF2s + 2 
index_barb = nsamples - 2 + 9
index_neom = nsamples - 1 + 9

stats.write('contig\t'         # chromosome name
            'position\t'       # bp position
            'indivs\t'         # individuals with data at this snp
            'med_reads\t'	   # median read depth per individual
            'RR\t'             # count of individuals with genotype RR
            'RA\t'             # count of individuals with genotype RA
            'AA\t'             # count of individuals with genotype AA
            'q\t'              # allele frequency of Ref allele
            'HWE_chisq\t'      # test for a deviation from hwe freqs
            'HWE_pval\t'       # significance for a deviation from hwe freqs
            'bar_geno\t'      # genotype call for hirsutus at this snp
            'neo_geno\n')     # genotype call for smallii at this snp

for line in in1:
    bar_geno = '0'
    neo_geno = '0'
    cols = line.replace('\n', '').split('\t')

    # Pass over metadata and header lines
    if len(cols) < 2:
        pass
    elif cols[0] == "#CHROM":
        pass

    else:
        tig = cols[0]
        pos = cols[1]
        totREF = 0
        totALT = 0
        calls = 0
        reads = []
        RRcount = 0
        RAcount = 0
        AAcount = 0

        for j in range(9, 9+nsamples):
            if cols[j] != './.':
                calls += 1
                fields = cols[j].split(':')  # e.g., 0/0:0,51,176:17:17,0:43

                if fields[0] != './.':
                    alleleDepths = fields[1].split(',')
                    refDepth = int(alleleDepths[0])
                    altDepth = int(alleleDepths[1])
                    reads.append(float(refDepth + altDepth))

                    phreds = fields[4].split(',')
                    if phreds[0] == '0':
                        RRcount += 1
                        if j==index_barb:
                            bar_geno = 'R'
                        elif j==index_neom:
                            neo_geno = 'R'
                    elif phreds[1] == '0':
                        RAcount += 1
                        if j==index_barb:
                            bar_geno = 'H'
                        elif j==index_neom:
                            neo_geno = 'H'
                    elif phreds[2] == '0':
                        AAcount += 1
                        if j==index_barb:
                            bar_geno = 'A'
                        elif j==index_neom:
                            neo_geno = 'A'

        medianRPI = median(reads)
        freqRR = float(RRcount) / calls
        freqRA = float(RAcount) / calls
        freqAA = float(AAcount) / calls
        freqR = freqRR + 0.5 * freqRA
        freqA = freqAA + 0.5 * freqRA
        expRR = calls * freqR * freqR
        expRA = calls * 2 * freqA * freqR
        expAA = calls * freqA * freqA

        if expRR == 0:
            termRR = 0
        else:
            termRR = (((RRcount - expRR) ** 2) / expRR)

        if expRA == 0:
            termRA = 0
        else:
            termRA = (((RAcount - expRA) ** 2) / expRA)

        if expAA == 0:
            termAA = 0
        else:
            termAA = (((AAcount - expAA) ** 2) / expAA)

        chisq = termRR + termRA + termAA

        pvalue = 1.0 - chi2.cdf(chisq, 1)

        stats.write(tig + '\t' +
                    pos + '\t' +
                    str(calls) + '\t' +
                    str(medianRPI) + '\t' +
                    str(RRcount) + '\t' +
                    str(RAcount) + '\t' +
                    str(AAcount) + '\t' +
                    str(freqR) + '\t' +
                    str(chisq) + '\t' + 
                    str(pvalue) + '\t' +
                    bar_geno + '\t' + 
                    neo_geno + '\n')

stats.close()






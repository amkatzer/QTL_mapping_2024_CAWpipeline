# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 15:04:53 2017

@author: carrie

getting basic info from vcf files

Assumes a filtered vcf!
"""

from numpy import median
from numpy import mean
from numpy import sum
from scipy.stats import chi2

in1 = open('allPheno_100individ_Q30.best.vcf', 'rU')
outfile = open('allPheno_100individ.best.hwe.vcf', 'w')

nF2s = 386
nsamples = nF2s + 2 
minq = 0.3
maxq = 0.7
minHWE = 0.0001

for line in in1:
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
            if cols[j].split(':')[0] != './.':
                calls += 1
                fields = cols[j].split(':')  
                
                alleleDepths = fields[1].split(',')
                refDepth = int(alleleDepths[0])
                altDepth = int(alleleDepths[1])
                reads.append(float(refDepth + altDepth))
                
                phreds = fields[4].split(',') ####0/1:145,0,80:23:6,17 example from bcftools output ###genotypecalls:phreds;readdepth:allelefreq

                if phreds[0] == '0':
                    RRcount += 1
                elif phreds[1] == '0':
                    RAcount += 1
                elif phreds[2] == '0':
                    AAcount += 1

        meanRPI = mean(reads)
        sumRPI = sum(reads)
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

        if pvalue >= minHWE and (minq <= freqR <= maxq):
            outfile.write(line)

outfile.close()

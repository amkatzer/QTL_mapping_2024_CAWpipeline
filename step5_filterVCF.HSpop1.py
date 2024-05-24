# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 18:04:58 2017

@author: carrie

this assumes the final two samples in the vcf are the parents.
also assumes already filtered for MQ > 30
"""
invcf = open('all_phenoed_NB_Q30.vcf', 'rU')
outvcf = open('20240523_CAWfilt_1individ_NB.vcf', 'w')
outinfo = open('filtreads_Q30_1individ.NB.info.txt', 'w')

outinfo.write('chr\tbp\tcount_F2s\n')

nF2s = 386 # number of F2s in the vcf file
minF2 = 1 # minimum number of F2s that need to have a particular snp to include it in the dataset
nsamples = nF2s + 2 # add 2 for the parents 


###################################################################

for line in invcf:
    cols = line.replace('\n','').split('\t')  # we are saving each line into a list

    # skipping header
    if len(cols) < 2:
        outvcf.write(line)
    elif cols[0] == '#CHROM':
        outvcf.write(line)

    else:
        if len(cols[4]) > 1: # if there is greater than 2 alt alleles, we skip because not biallelic
            pass
        else:

            # first check for fixed differences
            hirs = cols[nsamples-2 + 9] # pulling data from the hirsutus column
            smal = cols[nsamples-1 + 9] # pulling data from the smallii column

            if hirs.split(':')[0] != './.' and smal.split(':')[0] != './.': # if hirs and smal each have data
                if hirs.split(':')[0] != '0/1' and smal.split(':')[0] != '0/1': # if neither are hets
                    if hirs.split(':')[0] != smal.split(':')[0]: # if they are not the same type of homozygote

                        # seeing how many F2s have data
                        count_f2s = 0 # counter for the number of F2s with data
                        for j in range(9, 9+ nF2s):
                            if cols[j].split(':')[0] != './.':
                                count_f2s += 1

                        # is there data for at least 46 F2s?
                        if count_f2s >= minF2:
                            outvcf.write(line)
                            outinfo.write(cols[0] + '\t' + cols[1] + '\t' + str(count_f2s) + '\n')
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass

outvcf.close()
outinfo.close()


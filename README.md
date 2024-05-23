# QTL_mapping_2024_CAWpipeline
QTL mapping of Penstemon barbatus x Penstemon neomexicanus F2 population using Carrie Wessinger's protocol starting after the concatenation steps

### Step 1: concatenate files together
Concatenate any sequencing files together from multiple plates. Each sample should have the format "NB#" without a zero before the number (e.g. NB001 is incorrect, NB1 is correct).

Script to run for F2s: cat.fastq.sh

Code to run for parents: (bash, add a 'z' to the start of the species name to induce the parents to be alphabetically at the bottom of all of the samples including F2s)

```
find . -name "bar*" -print0 | xargs -0 cat -- > ../fastq_cat_output/zbarbatus.cat.fastq.gz

find . -name "neo*" -print0 | xargs -0 cat -- > ../fastq_cat_output/zneomexicanus.cat.fastq.gz
```

### Step 2: remove any samples from concatenation script output where there is no sequencing data. 
Delete anything with no data in the file to clean up the folder. I used the winscp gui to make this easy by sorting by sample size

### Step 3: make a sample list from the remaining outputs from the concatenation script
I used winscp gui again and copied the names (highlight all files -> right click -> file names -> copy to clipboard)

Save into a text file called testNames.txt.

### Step 4: use the python script to make scripts to run bwa-mem, samtools, and picard (adds read group info) for each sample
MAKE SURE YOU DO NOT GET EMAILS FOR ALL OF THESE SCRIPTS! They all run at the same time and I don't want to delete 600 emails. Can always check to see if they are running with squeue -u a681k477

```
mkdir bwa_shell_scripts
cd bwa_shell_scripts
module load python
python bwa_picard_samtools.py

#make directories
mkdir ../bwamem_output
mkdir ../picard_ARRG
mkdir ../picard_ARRG/stats

sh cg.cep.sh #this runs all the scripts
```

Check for errors in the log files:
grep "error" *.log

### Step 5: Run GATK UnifiedGenotyper

You will need to make a conda environment for gatk/3.8.0.

Edit script to add in your samples. Script: gatk_NB.sh

-> an easy way to make the -I NB#.bam\ for every sample is copying the sample names in winscp as described above into excel. In a column to the left of the sample names, write " -I " (be sure to include the spaces, but not the "). In the column to the right, write a backslash. Then in another column use the function concat() with the three columns. Copy and paste to get values rather than the function in a new column. 

Argument for memory is a Java argument. It looks like this in the script: -Xmx15G which indicates 15 GB mem

Check number of snps in file: grep -v "#" combined.vcf | wc -l

Get samples names from vcf into a text file: bcftools query -l input.vcf > samples.txt

### Step 5 extra! Checking the format of the output vcf

For all of the python scripts downstream, check what the format is of your initial VCF output file. If you do head -n 50 NB.vcf you should see something like this: GT:AD:DP:GQ:PL. You will have to change the called order of these in all the downstream scripts if yours is in a different order. 

Filter on Quality score
```
vcftools --vcf NB.vcf --out NB_Q30 --minQ 30 --stdout > NB_Q30.vcf
```

## Steps 6-10 will require you to edit the filenames & number of F2s within each of the scripts

### Step 6: Subset vcf to those that are phenotyped

Make a text file of all the phenotyped F2s + the parents that have the same ID as the vcf header.

bcftools view -S all_phenoed.txt -o all_phenoed.vcf NB.vcf

You may get an error such as: Error: subset called for sample that does not exist in header: "../bwamem_output/NB711.bam". Use "--force-samples" to ignore this error.

Copy the list of all genotyped IDs in the header of the vcf into the first column of an excel file.

Copy the list of all the phenotyped IDs into the second column of an excel file.

Run this in the third columm to see what has been included in the pheno that haven't been genotyped: =IF(COUNTIF(A:A,B2) > 0, TRUE, FALSE) 
  - Code from Dan Smith, excel wizard

Delete the FALSE values from the text file of all the phenotyped F2s + the parents (keep the parents even if they haven't been phenotyped)

run again with the new list: bcftools view -S all_phenoed.txt -o all_phenoed.vcf NB.vcf

### Step 7: Filtering out snps without substitutions between the two parents or with heterogeneity within a parent; filter by number of F2s who have to contain the allele

get number of F2s: wc -l samples.txt

subtract 2 from the output and edit the python script

python step5_filterVCF.HSpop1.py

### Step 8: find best snp then thin

Script to run first is finding the best snp: find.bestsnp.py

Second script to run: thin.bestsnp.py

### Step 9: Calculate HWE

Calculate the values: hwe.calc.py

Graph the allele frequency
```
module load R
R
data <- read.table("allPheno_100individ.best.vcf.info.txt", header=T)
pdf("20240522_allele_freq.pdf")
hist(data$q)
dev.off()
quit()
```

filter by HWE & allele frequency: filtervcf_q_hwe.py

### Step 10: make files for lepmap

Edit your samples file in a new file that is only the F2s.

Run script: 20240513_generate_lepmap.AMK.py

### Step 11: Lepmap

You will need to put a folder of the lepmap software on the cluster where you can call it. Mine is in scratch currently.

All of the steps require java, so run module load java.

Code from Trinity Depatie!

```
#put lepmap into a folder on the cluster (probably scratch)

module load java

###Step 1: calling genotypes in parents 
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin ParentCall2 data=./allPheno_100individ.best.lepmap.input.txt removeNonInformative=1 > filter.parentcall.txt

###Step 2: filtering markers generated in previous command 
###dataTolerance paramter affects segregation distortion and missingLimit affects excess number of missing genotypes
###single family data should have a smaller value for data tolerance
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin Filtering2 data=./filter.parentcall.txt dataTolerance=0.000001 > filter.filteredcalls.txt

###Step 3: save names and positions of contigs to a file
cat filter.filteredcalls.txt|cut -f 1,2|awk '(NR>=7)' > filter.snps.txt

###Once you have the snps file, you want to download that off the cluster and open in excel. Add a column at the end which is named "snp_number" or something and start with 1 and allow excel to fill in the rest of the ordered numbers until the end of the document. 
###Now that you have this additional column in the excel, you are going to make a new excel document that you will paste into. From the original document, copy the numbers in the "snp_number" column that correspond with each LG. Paste these numbers into a new excel document (different document for each LG).
###Now you should have 8 different files (in the below command, you can see that I named mine "chr<_>.manual.order.126F2.txt"). The file corresponding with the first linkage group should only contain 1 column with row values starting at 1 and continuing until the # of snps you have on that LG. 

###AFTER MANUAL REORDERING### 
###To manually order markers, download the SNPs file previously generated. Then add a column (snp_num) that provides a number that corresponds with each marker starting with 1. Filter these files to only contain the snp num that corresponds with the marker. These new txt files should be saved as chr*.manual.order.<additional_info>.txt in order to use my script with ease.
###Place markers in proper physical position on each chromosome 
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr1.manual.order.txt data=filter.filteredcalls.txt chromosome=1 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr1.manually.ordered.txt
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr2.manual.order.txt data=filter.filteredcalls.txt chromosome=2 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr2.manually.ordered.txt
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr3.manual.order.txt data=filter.filteredcalls.txt chromosome=3 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr3.manually.ordered.txt
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr4.manual.order.txt data=filter.filteredcalls.txt chromosome=4 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr4.manually.ordered.txt
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr5.manual.order.txt data=filter.filteredcalls.txt chromosome=5 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr5.manually.ordered.txt
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr6.manual.order.txt data=filter.filteredcalls.txt chromosome=6 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr6.manually.ordered.txt
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr7.manual.order.txt data=filter.filteredcalls.txt chromosome=7 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr7.manually.ordered.txt
java -cp /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin OrderMarkers2 evaluateOrder=chr8.manual.order.txt data=filter.filteredcalls.txt chromosome=8 grandparentPhase=1 outputPhasedData=1 sexAveraged=1 improveOrder=0 > chr8.manually.ordered.txt

###Convert phased marker format into a useful format with genotype data for each chromosome with MANUALLY ORDERED DATA
###The map2genptypes.awk script needs to be downloaded from a third party website. This script did not come with my version of Lepmap3.
awk -vchr=1 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr1.manually.ordered.txt > chr1.manually.ordered.with.genos.txt
awk -vchr=2 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr2.manually.ordered.txt > chr2.manually.ordered.with.genos.txt
awk -vchr=3 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr3.manually.ordered.txt > chr3.manually.ordered.with.genos.txt
awk -vchr=4 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr4.manually.ordered.txt > chr4.manually.ordered.with.genos.txt
awk -vchr=5 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr5.manually.ordered.txt > chr5.manually.ordered.with.genos.txt
awk -vchr=6 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr6.manually.ordered.txt > chr6.manually.ordered.with.genos.txt
awk -vchr=7 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr7.manually.ordered.txt > chr7.manually.ordered.with.genos.txt
awk -vchr=8 -f /kuhpc/scratch/hileman/a681k477/QTL_working/20240306_vcftools/lepmap_working/lepmap3_software/bin/map2genotypes.awk chr8.manually.ordered.txt > chr8.manually.ordered.with.genos.txt

###Convert Lepmap3 marker numbers back to genomic coordinates
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr1.manually.ordered.with.genos.txt > chr1.manually.ordered.with.genos.mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr2.manually.ordered.with.genos.txt > chr2.manually.ordered.with.genos.mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr3.manually.ordered.with.genos.txt > chr3.manually.ordered.with.genos.mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr4.manually.ordered.with.genos.txt > chr4.manually.ordered.with.genos.mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr5.manually.ordered.with.genos.txt > chr5.manually.ordered.with.genos.mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr6.manually.ordered.with.genos.txt > chr6.manually.ordered.with.genos.mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr7.manually.ordered.with.genos.txt > chr7.manually.ordered.with.genos.mapped.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' filter.snps.txt chr8.manually.ordered.with.genos.txt > chr8.manually.ordered.with.genos.mapped.txt
```

Graph the linkage map:
```
module load R
R

chr1 <- read.table("./chr1.manually.ordered.with.genos.mapped.txt")
chr2 <- read.table("./chr2.manually.ordered.with.genos.mapped.txt")
chr3 <- read.table("./chr3.manually.ordered.with.genos.mapped.txt")
chr4 <- read.table("./chr4.manually.ordered.with.genos.mapped.txt")
chr5 <- read.table("./chr5.manually.ordered.with.genos.mapped.txt")
chr6 <- read.table("./chr6.manually.ordered.with.genos.mapped.txt")
chr7 <- read.table("./chr7.manually.ordered.with.genos.mapped.txt")
chr8 <- read.table("./chr8.manually.ordered.with.genos.mapped.txt")

pdf("physical_linkage_plot_all_chr.pdf")
par(mfrow = c(3,3))
plot(chr1$V2, chr1$V4, xlab="bp", ylab="cM", main="Chr1")
plot(chr2$V2, chr2$V4, xlab="bp", ylab="cM", main="Chr2")
plot(chr3$V2, chr3$V4, xlab="bp", ylab="cM", main="Chr3")
plot(chr4$V2, chr4$V4, xlab="bp", ylab="cM", main="Chr4")
plot(chr5$V2, chr5$V4, xlab="bp", ylab="cM", main="Chr5")
plot(chr6$V2, chr6$V4, xlab="bp", ylab="cM", main="Chr6")
plot(chr7$V2, chr7$V4, xlab="bp", ylab="cM", main="Chr7")
plot(chr8$V2, chr8$V4, xlab="bp", ylab="cM", main="Chr8")
dev.off()
quit()
```

### Step 12: Generate input for r/QTL

script to run: generate_rqtl_input_file.py

### Step 13: r/QTL

Code from Haylee Nedblake!


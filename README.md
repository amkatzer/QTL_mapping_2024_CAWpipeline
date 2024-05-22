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

argument for memory is a Java argument. It looks like this in the script: -Xmx15G which indicates 15 GB mem




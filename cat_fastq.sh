#!/bin/bash
#SBATCH --job-name=cat_fastq
#SBATCH --partition=sixhour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a681k477@ku.edu
#SBATCH --time=06:00:00
#SBATCH --output=cat_fastq_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb

mkdir ../fastq_cat_output/

for i in $(cat F2_list.txt)
do
F2name="NB"$i
current_F2_pathway=$(find . -name "$F2name.*")
echo $F2name $current_F2_pathway
find . -name $F2name.* -print0 | xargs -0 cat -- > ../fastq_cat_output/$F2name.cat.fastq.gz
done

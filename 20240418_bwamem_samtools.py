commandsFile = open("cg.cep.sh", 'w')
sampleNames = open("testNames.txt", "r")
data = sampleNames.readlines()

for rp in data:
	rp = rp.strip()
	new_rp = rp.replace(".cat.fastq.gz", "")
	sampleName = "t_"+str(new_rp) #

	shellScriptName = 'r%s.sh' % (sampleName)
	shellScript = open(shellScriptName, 'w' )

	commandsFile.write('sbatch %s \n' % (shellScriptName))

	shellScript.write("#!/bin/bash\n" )
	shellScript.write("#SBATCH --job-name=%s\n" % (sampleName))
	shellScript.write("#SBATCH --partition=sixhour       \n" )
	shellScript.write("#SBATCH --time=0-5:59:00        \n" )
	shellScript.write("#SBATCH --mem=5gb           \n" )
	shellScript.write("#SBATCH --mail-type=NONE          \n" )
	shellScript.write("#SBATCH --mail-user=a681k477@ku.edu    \n" )
	shellScript.write("#SBATCH --ntasks=1                   \n" )
	shellScript.write("#SBATCH --cpus-per-task=5            \n" )
	shellScript.write("#SBATCH --output=k_%s.log\n\n\n" % sampleName)

	shellScript.write("module load bwa \n" )
	shellScript.write("mem -t 5 Pbar.2022.LG.fa ../../../fastq_cat_output/"+sampleName+".cat.fastq.gz | samtools sort -o ../bwamem_output/"+sampleName+".bam \n" )
	shellScript.write("picard AddOrReplaceReadGroups I=../bwamem_output/"+sampleName+".bam O=../picard_ARRG/"+sampleName+".bam SO=coordinate RGID="+sampleName+" RGLB=../bwamem_output/"+sampleName+".bam RGPL=illumina RGPU=../bwamem_output/"+sampleName+".bam RGSM=../bwamem_output/"+sampleName+".bam VALIDATION_STRINGENCY=LENIENT \n")
	shellScript.write("samtools index ../picard_ARRG/"+sampleName+".bam ../picard_ARRG/"+sampleName+".bai \n")
	shellScript.write("samtools stats ../picard_ARRG/"+sampleName+".bam > picard_ARRG/stats/"+sampleName+".bam.stats")

sampleNames.close()
commandsFile.close()

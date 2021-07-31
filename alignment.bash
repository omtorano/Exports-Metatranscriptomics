#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-03:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=omtorano@email.unc.edu
#SBATCH -J salmon_trin
#SBATCH -o salmon_trin.%A.out
#SBATCH -e salmon_trin.%A.err

#ran in Job Wall-clock time: 1-13:31:16
#Memory Utilized: 14.01 GB
#to create index

module load salmon

indir=/proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/Reads/pairedendtrimmed

input1=`ls ${indir}/*-1*R1*gz | awk -F 'R1' '{print $1}'`
input4=`ls ${indir}/*-4*R1*gz | awk -F 'R1' '{print $1}'`


#First need to create index, then map
salmon index -i /pine/scr/o/m/omtorano/alignment/trinity/assemblyindex \--transcripts /pine/scr/o/m/omtorano/annotation/trinity/fastanno_out/clustered_assembly_depth1_trinity.fasta -k 31
for s in `echo $input1`; do 

salmon quant -l A -i /pine/scr/o/m/omtorano/alignment/trinity/assemblyindex \
	-1 ${input1}R1*.gz \
	-2 $indir/R2*.gz \
	-p 24 \
	-o /pine/scr/o/m/omtorano/alignment/trinity/${input1}_quants_trinity
done

salmon index -i /pine/scr/o/m/omtorano/alignment/trinity/assemblyindex \--transcripts /pine/scr/o/m/omtorano/annotation/trinity/fastanno_out/clustered_assembly_depth4_trinity.fasta -k 31
for s in `echo $input4`; do 

salmon quant -l A -i /pine/scr/o/m/omtorano/alignment/trinity/assemblyindex \
	-1 ${input1}R1*.gz \
	-2 $indir/R2*.gz \
	-p 24 \
	-o /pine/scr/o/m/omtorano/alignment/trinity/${input1}_quants_trinity
done
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=00-4:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=5
#SBATCH --array=1-21
#SBATCH -J trim
#SBATCH -o trim.%A_%a.out
#SBATCH -e trim.%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=omtorano@email.unc.edu

#load necessary modules for trim_galore, auto loads python, cutadapt
module load trim_galore
module load pigz


#set in directory to where raw reads are
indir=/pine/scr/o/m/omtorano/exports/reads
outdir=/pine/scr/o/m/omtorano/exports/trimmed_reads

#create out directory
echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi

RUN=${SLURM_ARRAY_TASK_ID}
#cuts path file at R1, prints everything before
input=`ls ${indir}/*R1*.gz | awk -F 'R1' '{print $1}'| sed -n ${RUN}p`


echo -e "\nRun ID: ${RUN}"
echo -e "\nSample: ${input}"


trim_galore -j 4 \
	--paired ${input}R1* ${input}R2* \
	-o  ${outdir}
	




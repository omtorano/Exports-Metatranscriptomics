#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=02-0:00:00
#SBATCH --mem=60G
#SBATCH --ntasks=42
#SBATCH -J qc
#SBATCH -o qc.%A.out
#SBATCH -e qc.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=omtorano@email.unc.edu

module load fastqc

#set raw_reads directory to where raw reads are
raw_reads=/pine/scr/o/m/omtorano/exports/reads
#set trim_reads directory to be where trimmed reads are
trim_reads=/pine/scr/o/m/omtorano/exports/trimmed_reads
#set out directory to where fastqc output should go, can make 2 directories to keep raw and trimmed output seperate
outdir=/pine/scr/o/m/omtorano/exports/fastqc_output

#create ourdirectory
echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi

#run fastqc on untrimmed raw reads, reads will run simultaneously, set threads to # of reads
fastqc -t 42 $raw_reads/* -o ${outdir}
#run fastqc on trimmed raw reads, using *.gz searches for anything that ends with .gz, excludes trimming reports
fastqc -t 42 $trim_reads/*.gz -o ${outdir}

#create multiqc file for all fastqc output
module load multiqc
cd ${outdir}
multiqc . 
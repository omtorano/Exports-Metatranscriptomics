#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 3-00:00:00
#SBATCH -J cdhit_trinity
#SBATCH -o cdhit_trinity.%j.out
#SBATCH -e cdhit_trinity.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=omtorano@email.unc.edu
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1

indir=/proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/assemblies/
outdir=/pine/scr/o/m/omtorano/cdhit/
# ifn is the combined (concatenated) assemblies
input1=`ls ${indir}/*-1*fasta`
input4=`ls ${indir}/*-4*fasta`
ifn1="${outdir}/mega_assembly_depth1_trinity_renamed.fasta"
ifn4="${outdir}/mega_assembly_depth4_trinity_renamed.fasta"
ofn1="${outdir}/clustered_assembly_depth1.fasta"
ofn4="${outdir}/clustered_assembly_depth4.fasta"

echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi

echo "Checking if ${ifn1} exists ..."
if [ ! -f "${ifn1}" ]
then
    echo "Create combined assembly ... ${ifn1}"
    cat ${input1} > ${ifn1}
else
    echo " ... exists"
fi

echo "Checking if ${ifn4} exists ..."
if [ ! -f "${ifn4}" ]
then
    echo "Create combined assembly ... ${ifn4}"
    cat ${input4} > ${ifn4}
else
    echo " ... exists"
fi



# load default cdit
module load cdhit

echo "${ifn1}"
echo "${ofn1}"

cd-hit-est \
 -i "${ifn}" \
 -o "${ofn}" \
 -c .98 -n 10 -d 100 \
 -T ${SLURM_CPUS_PER_TASK} \
 -M 40000 \

echo "${ifn4}"
echo "${ofn4}"

cd-hit-est \
 -i "${ifn4}" \
 -o "${ofn4}" \
 -c .98 -n 10 -d 100 \
 -T ${SLURM_CPUS_PER_TASK} \
 -M 40000 \

# --------------------- 
# sacct -j $SLURM_JOB_ID --format='JobID,user,elapsed, cputime, totalCPU,MaxRSS,MaxVMSize,ncpus,NTasks,ExitCode'

scontrol show job $SLURM_JOB_ID

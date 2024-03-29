#!/bin/bash
#SBATCH -N 1
#SBATCH -t 6-00:00:00
#SBATCH -J trinity
#SBATCH -o trinity.%j_%a.out
#SBATCH -e trinity.%j_%a.err
#SBATCH --array=1-21%3
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=omtorano@unc.edu
#SBATCH --mem=500G
#SBATCH --cpus-per-task=3
#SBATCH --ntasks=1


module load trinity
#set in directory to where trimmed reads are stored
indir=/pine/scr/o/m/omtorano/exports/trimmed_reads
RUN=${SLURM_ARRAY_TASK_ID}
input=`ls ${indir}/*R1*gz | awk -F 'R1' '{print $1}'| sed -n ${RUN}p`
a=`basename ${input}`
#make unique out directory for each sample by amending sample name to trinity outdir !!out directory must contain 'trinity'!!
outdir=/pine/scr/o/m/omtorano/exports/assembly/${a}trinity
#hpc_cmds_GridRunner.pl is stored here
griddir=/proj/marchlab/tools/HpcGridRunner

echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi

 echo "Checking completion ..."
if [ ! -f ${outdir}.Trinity.fasta ]
then
echo "creating ${outdir}.Trinity.fasta"
Trinity \
	--seqType fq \
	--max_memory 500G \
	--left ${input}R1*gz \
	--right ${input}R2*gz \
	--CPU ${SLURM_CPUS_PER_TASK} \
	--grid_exec "${griddir}/hpc_cmds_GridRunner.pl --grid_conf ${griddir}/hpc_conf/trinity_pipeline_longleaf_grid.conf -c" \
	--full_cleanup \
	--NO_SEQTK \
	--output ${outdir}
fi

 echo "Checking completion ..."
if [ ! -f ${outdir}.Trinity.fasta ]
then
echo "${outdir}.Trinity.fasta not complete... resubmitting"

cat<<EOF>${a}trinity_rerun.bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH -J ${a}trinity_rerun
#SBATCH -o ${a}trinity_rerun.%j.out
#SBATCH -e ${a}trinity_rerun.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

module load trinity

Trinity \
	--seqType fq \
	--max_memory 30G \
	--left ${input}R1*gz \
	--right ${input}R2*gz \
	--CPU ${SLURM_CPUS_PER_TASK} \
	--grid_exec "${griddir}/hpc_cmds_GridRunner.pl --grid_conf ${griddir}/hpc_conf/trinity_pipeline_longleaf_grid.conf -c" \
	--full_cleanup \
	--NO_SEQTK \
	--output ${outdir}
	/nas/longleaf/apps/trinity/2.8.6/trinityrnaseq-2.8.6/util/TrinityStats.pl ${outdir}/Trinity.fasta
 
echo "Checking completion ..."
if [ ! -f ${outdir}.Trinity.fasta ]
then
sbatch ${a}trinity_rerun.bash
else
    echo " ${outdir}.Trinity.fasta complete... Delete Trinity directory."
	/nas/longleaf/apps/trinity/2.8.6/trinityrnaseq-2.8.6/util/TrinityStats.pl ${outdir}.Trinity.fasta > ${outdir}.Trinity.stats
	rm -r ${outdir}
fi


EOF

sbatch ${a}trinity_rerun.bash
else
    echo " ${outdir}.Trinity.fasta complete... Delete Trinity directory."
	/nas/longleaf/apps/trinity/2.8.6/trinityrnaseq-2.8.6/util/TrinityStats.pl ${outdir}.Trinity.fasta > ${outdir}.Trinity.stats
	rm -r ${outdir}
fi

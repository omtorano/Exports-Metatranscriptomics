#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=omtorano@email.unc.edu
#SBATCH -J diamond_trinity
#SBATCH -o diamond_trinity.%j.out
#SBATCH -e diamond_trinity.%j.err

module load diamond


indir=/proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/cdhit/

outdir=/pine/scr/o/m/omtorano/annotation/

#diamond makedb --in /proj/marchlab/data/phylodb/phylodb_1.076.pep.fa -d phylodb
#diamond makedb --in /nas/longleaf/data/KEGG/KEGG/genes/fasta/genes.pep.fasta -d keggdb

samples='clustered_assembly_depth1.fasta clustered_assembly_depth4.fasta'

for s in `echo $samples`; do            

#Need to blast against phyloDB and KEGG

diamond makedb --in /proj/marchlab/data/phylodb/phylodb_1.076.pep.fa -d phylodb

diamond blastx -d /proj/marchlab/data/phylodb/diamond_db/phylodb \
	-q $indir/${s} \
	-o $outdir/${s}phyloDB.m8 \
	-p 12 -e 0.000001 -k 1
	#-p = threads, -e =Maximum expected value to report an alignment (default=0.001), -k = max target seqs number (how many seqs/ puery to report alignments)

diamond makedb --in /nas/longleaf/data/KEGG/KEGG/genes/fasta/genes.pep.fasta -d keggdb

diamond blastx -d keggdb \
	-q $indir/${s} \
	-o $outdir/${s}kegg.m8  \
	-p 12 -e 0.000001 -k 1
	
	
done
# Marchetti-lab-metatranscriptomics
Annotated metatranscriptomics pipeline
adapted from https://github.com/Lswhiteh/diatom-metatranscriptomics
# Overview
This pipeline was used for the 2020 High Yield EXPORTS metatranscriptomics RNA samples, but every effort has been made to make these scripts easy to understand and modify for future projects. 

The general steps this pipeline follows are:
  - Trimming
  - Quality control
  - Assembly
  - Annotation 
  - Alignment
  - Differential expression

## Before Getting started
Once RNA extraction is complete, libraries have been prepped and sequenced, and sequences have been delivered the following steps are necessary to access and analyze seqs

- Get on longleaf (UNC's high performance computing platform)
	- Follow the instructions here https://its.unc.edu/research-computing/request-a-cluster-account/
		- You'll need the following info, -> are my suggestions
			- Your name and onyen: 
			- Your Department: -> Marine Sciences or MASC
			- Phone number to reach you while your jobs are running (if we need to): 
			- Your Email address:
			- Your preferred shell (bash or tcsh): -> bash
			- Faculty sponsor name: -> Adrian Marchetti
			- Faculty sponsor onyen: -> amarchet
			- Type of subscription: Longleaf or dogwood -> longleaf
			- A description of the type of work you will do on the cluster -> bioinformatics/RNA seq/etc. something like that
	- After you get longleaf access you may need to email ITS to get access to the Marchetti lab directory, in your email request access to /proj/marchlab and cc Adrian
- Find a platform for accessing longleaf that works for you https://its.unc.edu/research-computing/getting-logged-on/
	- There are many ways to do this. Some options include using longleaf desktop https://ondemand.rc.unc.edu/pun/sys/dashboard/ (Sarah's preferred method), doing everything (navigating directories, editing files, etc.) from a terminal (Johnson's preferred method), doing hybrid terminal and file manager (my preferred method, on Windows I use a GitBash and WinSCP file manager, it works for me). 
- Using longleaf off campus requires a VPN https://help.unc.edu/sp?id=kb_article_view&sysparm_article=KB0010155&sys_kb_id=719db1eddb3fa41070551ffa689619eb
- Make a new directory within /proj/marchlab/projects for your project, or new subdirectory within one of the existing project directories. Then download sequences if sequencing was done by an outside service (not UNC's High Throughput Sequencing facility HTSF)
	- For example - when I got sequences for the High Yield 2020 EXPORTS project I made a directory within the existing Exports project folder called /HighYield2020. The sequencing was completed at Genewiz, who sent email detailing how to transfer files with sftp. This is an example of what this would look like (JUST AN EXAMPLE, DO NOT RUN):
```
### log into longleaf with your email

ssh omtorano@longleaf.unc.edu

### password prompt will come up automatically, you will not be able to see characters as you type them

Password: <enter password here>

### cd changes directory

cd /proj/marchlab/projects/EXPORTS/metatranscriptomics/

### mkdir makes a directory, here I am making the directory HighYield2020 and another directory within that called Reads on the following line

mkdir HighYield2020
mkdir HighYield2020/Reads

### the beginning of the directions from Genewiz - connect to their server via sftp using the email in the instructions

sftp a_marchetti_lab_gmail@sftp.genewiz.com

### password prompt will come up automatically, enter the password they sent via email

password: <enter password genewiz sent here>

### lcd sets the local directory (where you want your reads to go), does not change working directory (does not change where you are)

lcd /proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/Reads

### ls lists items in the current working directory, in this case I just logged onto Genewiz so I am seeing the list of directories
### I have access to - likely there will be one and it will have the name of the Genewiz project number

ls

### change directory into whatever directory was just listed by ls

cd <genewiz project folder name that just came up from ls command>

### mget transfers files from working directory to whatever you set as your local directory. The '*' is a Linux wildcard that here 
### means 'get all the stuff'

mget *

### end the sftp connection to transfer yourself back to longleaf

quit #Is this the command that ends sftp connection?
```
## Getting Started
Rules for the Marchetti Lab /proj space
- Do not work in the proj space, it is only for storage of 'final product' files from each stage in pipeline. We have a history of running out of space, storage availability can be tracked here https://rc-storage-info.its.unc.edu:32000.
- Instead work in scratch (/pine/scr/o/n/onyen, for example my path is /pine/scr/o/m/omtorano) or home directory (/nas/longleaf/home/onyen). Scratch has tons of space but inactive files will be deleted by ITS after ~30days, they send an email before deletion.
- Use the rsync command to move files, mv deletes files from their original location and if interrupted will result in data loss.
- Use the rm (remove) command CAREFULLY, permanently deletes files.
- The recommended directory organization is to have directory for the project and separate sub directories for code, reads, assemblies, annotations, and alignment (optional trimmed reads and fastqc output).
- See final section of page for useful stuff.

# Trimming with trim_galore

Tl;dr
```
sbatch trim_galore.bash
```
The details...

For trimming I used the tool trim_galore https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
This tool removes low quality reads and auto detects and removes common adapters. Here I will explain the trimgalore.bash script, to download a version that can be run (with minor edits to file path names etc.) see files listed above or git clone. Copy and pasting individual sections of code will not work correctly, this must be run as a script.

1) Longleaf is a big computer with lots of nodes, when you log on to longleaf you, and everyone else, are on the login node. Anything you do on the login node takes up computing resources, if longleaf is being really slow it could be because there are tons of people using it. Due to the limited resources, anything more than really resource light commands (i.e. moving files) should be done as a 'job'. A 'job' is basically what you call the commands you have saved in a file called a script once you ask it to be run. The way you ask for longleaf to run a job is by using the 'sbatch' command followed by the name of your script. Sbatch is a command that communicates with SLURM, the job scheduler that is used on longleaf. SLURM is just one of a few scheduling tools used by high performance computing platforms, but it is common enough that there are a lot of resources online. SLURM basically manages the logistics of allocating compute time, memory, nodes etc. See https://its.unc.edu/research-computing/techdocs/getting-started-on-longleaf/#Job%20Submission. To tell SLURM the resources you need the following paragraph goes at the beginning of a script:
```
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
```
While '#' normally indicates lines within code that are not read, that isn’t the case here. The first line is called a shebang and for this script tells longleaf this is a bash script, other scripts (python, etc.) might look different. Below this are the SLURM options, here it's requesting a general node (instead of bigmem, snp, etc.), 1 node, 4 hours of time, 10 gigs of memory, 5 tasks, indicating that we are requesting an array and giving the length of the array (here I have 21 samples), assigning the job name to be trim, assigning the name of the '.out' file to be trim.jobnumber_arraynumber.out, same for the error file but .err, requesting longleaf send an email when the job starts, ends, or fails, and gives the email address to send the notification to. Most (all?) functions when run create a 'output stream' and 'error stream', basically any text output from a function working properly will be part of the 'output stream', and errors will be part of the 'error stream'. Instead of printing to the terminal, these streams are passed into files that will appear in your working directory. This job creates 21 .out and 21 .err files, hence it is important to add the unique job number to the name. 

2) load necessary modules for trim_galore, auto loads python & cutadapt
```
module load trim_galore
module load pigz
```
Longleaf has many software tools available that can be used by adding them to your working space with 'module load' ('module add' does the same thing). Other useful module management commands: module avail - lists allllll available modules (modify by doing 'module avail r*' etc. to search for modules starting with r etc.), 'module list' - shows all the modules you have loaded in your working space (loading a module in a job does not mean it will be loaded in your login session), 'module rm <module name>' - removes specific module, 'module purge' - removes all modules loaded in your space. Most longleaf modules are stored in longleaf here: /nas/longleaf/apps/

3) set 'in' directory to where raw reads are, if you need to transfer reads from /proj run rsync -r (recursive, needed if transferring all files in directory) /from/path /to/path
```
rsync -r /proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/Reads /pine/scr/o/m/omtorano
indir=/pine/scr/o/m/omtorano/exports/reads
outdir=/pine/scr/o/m/omtorano/exports/trimmed_reads
```
Transferring files is not necessary, a script run out of scratch can access files in any directory (assuming the user has permission). In the above case I am basically making a copy of my raw reads in my scratch space. This might be desirable if multiple users are using files and they are being moved, or you want to make a copy to ensure nothing happens to the original read files. 
4) This creates the out directory if it does not already exist. Echo is a Linux command which basically means print - i.e. the line "checking if out directory exists" will print in your .out file. This is not necessary for creating the directory, just helpful to tell you what it’s doing. Below this is a for loop and conditional saying if out directory does not exist make it, if it does exist print "... exists". The "!" means does not exist and the -d indicates what follows will be a directory, in this case if outdir does not exist. 
```
echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi
```
5) Set variable called 'RUN' to = slurm array task id (https://slurm.schedmd.com/job_array.html). This is a very cool tool for submitting a bunch of jobs at once, more explanation below. Naming really necessary but makes the following line a bit cleaner.
```
RUN=${SLURM_ARRAY_TASK_ID}
```
6) Parse and define input as the sample name that will be input to trim_galore. ls ${indir}/&ast;R1&ast; will list all files in in directory that have R1 anywhere in the name, the * wildcard is used here to mean find R1 embedded in any string of characters before or after. &ast;R1 would only look for R1 occurring at the end of any string of characters, R1&ast; would only look for R1 at the beginning. This is done to get each sample name listed one time (there will be an R1 and R2 in the indir assuming paired end reads). The '|' character 'pipes' commands, the output of ls is given to awk. awk here is being used to cut path file name at R1, -F is the field separator which here is set to R1. This results in names being split before and after the R1, '{print $1}' prints the first element of the split name, now we have each unique sample ID listed one time. sed here is being used to isolate the 'changed' input names. The slurm array task id will run through a list of all of the unique sample names, to input them one at a time sed -n followed by p prints only the file being processed by the slurm array. sed and awk are useful Linux tools but there are multiple ways to accomplish this string parsing to get the unique sample names.
```
input=`ls ${indir}/*R1* | awk -F 'R1' '{print $1}'| sed -n ${RUN}p`
```
7) Echo (print) the run id (here 1-21) and sample name purely for user ease of comprehension, this has no effect on trimming. The -e flag here allows echo to interpret the backslash escapes, here used to print on a new line "\n", other uses include "\t" if you want to echo something but have it be tab separated.
```
echo -e "\nRun ID: ${RUN}"
echo -e "\nSample: ${input}"
```
8) Now the actual trimming part, -j tells it how many cores to use, --paired tells it these are paired end reads, ${input} puts the unique sample id being cycled through by the slurm array, and we add back on the R1 and R2 endings, using the * wildcard lets us not have to type the entire file name. I think due to python versions it will actually only run one core as written here.
```
trim_galore -j 4 \
	--paired ${input}R1* ${input}R2* \
	-o  ${outdir}
```

9) To run this script download it or git clone, edit the path files and email address, put the script into your working directory (i.e. your scratch space /pine/scr/o/n/onyen) and type the following:
```
sbatch trim_galore.bash
```
If working on one project in scratch it would not be necessary to make further subdirectories for reads as is written here (reads is a subdirectory within exports in this example).

# Quality control with Fastqc
In this and following steps I will only go through new code elements that differ from trim_galore.
No array is used here, set tasks to # of samples. Due to the lack of array elements the a% part of .out and .err file names is unnecessary.
```
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

raw_reads=/pine/scr/o/m/omtorano/exports/reads        
trim_reads=/pine/scr/o/m/omtorano/exports/trimmed_reads
```
set out directory to where fastqc output should go, can make 2 directories to keep raw and trimmed output separate.
```
outdir=/pine/scr/o/m/omtorano/exports/fastqc_output

#create out directory
echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi

```
run fastqc on untrimmed raw reads and trimmed reads, reads will run simultaneously, set threads to # of reads.
```
fastqc -t 42 $raw_reads/* -o ${outdir}
#run fastqc on trimmed reads, using *.gz searches for anything that ends with .gz, excludes trimming reports
fastqc -t 42 $trim_reads/*.gz -o ${outdir}
```
create multiqc file for all fastqc output, '.' means 'current directory'.
```
cd ${outdir}
module load multiqc
multiqc . 
```
The multiqc step will result in an .html file that can be dragged and dropped into a browser. This file contains the compiled fastq reports from all trimmed and raw reads. If there are certain samples that have low quality or vary from the rest of the samples this may indicate that these sample need to be removed from further analysis. Full multiqc output for the High Yield 2020 samples is saved in HighYield2020_multiqc_report.html, below is an example of what it looks like... 
	
![image](https://user-images.githubusercontent.com/48129653/123166593-96a19000-d443-11eb-8f6b-e47c3c7fc0a1.png)


# Assembly with Trinity
The goal of assembly in this work flow is to create a de novo assembly from all samples in a 'group' which you can use to map back to in the alignment step. For Exports High Yield 2020, samples were collected at two timepoints from two different depths, 1 & 4. I am interested in comparing between timepoints, not depths, I am therefore creating an assembly of all depth 1 samples and all depth 4 samples. This way I can compare samples collected at a single depth, instead of between both depths. 

Trinity is one option for assembly tools, the other I will go through here (next section) is Spades. Trinity is a multipurpose tool that can do a lot more than assembly, read more here https://github.com/trinityrnaseq/trinityrnaseq/wiki. During assembly it runs through different phases (see wiki) that can be split into different jobs depending on how big an assembly you are trying to make. For this project the stats of my assemblies are below:

![image](https://user-images.githubusercontent.com/48129653/119878716-7809b100-bef8-11eb-97c4-d12eb5d70787.png)

Totals were calculated based on averages (total sequence length = average sequence length * number of samples), averages were calculated during multiqc step.

Assembly is where things start to get complicated, depending on how many sequences you have and the questions you are trying to answer there are many different tools and ways to do it within those tools - there are lots of studies comparing different tools. I tried MANY different methods for assembly with varying success. The script described here will assemble each individual sample. They will then be combined into a depth 1 & 4 assembly in the next step (cdhit). 
Some other things I tried: 
	- assembling samples by triplicate and assembling those in cdhit - this resulted in relativly low mapping rates and took a lot of time and memory. In fact, assembling triplicates required Trinity to be broken up into multiple scripts due to maxed out memory.
	- assembling by depth (all depth 1 samples in one script and all depth 4 in another), never finished withing alloted time and memory even when broken down by phase. 
	- using IDBA (never finished).
	- Spades individual assemblies (see next section)
	- Spades depth assemblies (see next section)
	
As of 6/23/3031 the newest version of Trinity available on longleaf is trinity/2.8.6 and the newest version of trinity available is trinity/2.12.0. There have been some requests for ITS to update trinity but it has not happened. Options for running trinity include:
	- using the available version on longleaf, easiest and what is described below 
	- using a docker or singularity image (https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#trinity_singularity), this should be easy but trinity does not currently support the use of the gridrunner tool (see below) on singularity
	- downloading and compiling the newest version of trinity onto longleaf. I did not spend a lot of time on this but it is not that easy and I did not get it to work. When trinity is added using 'module load' additional packages are automatically added and the user environment (view by typing 'env' into command line) is automatically modified so all the libraries etc. have the right paths. This would need to be altered manually if trinity is run from a manually compiled version.
	
Module loading trinity will automatically load bowtie/1.2.3, bowtie2/2.4.1, jellyfish/2.2.10, trinity/2.8.6, perl/5.18.2, samtools/1.12, and salmon/0.9.1. Longleaf auto loading other necessary packages is a common thing and good to know about for a few reasons. 1) It does not always auto load the correct versions of these software tools. When I first started using trinity an old version of samtools was autoloaded, I emailed ITS about this and they fixed it so the correct version loaded. 2) If you are using multiple tools there may be conflicting versions of auto loaded software, this is common with python and perl. An additional to the array flag is '%3', this will submit 3 jobs at once instead of submitting all 21. This is important because no directory, including scratch, has enough space to hold the intermediate files of more than ~6 trinity runs at once. I am going 3 here to be safe, if you reach your scratch space & file limit all jobs will fail.
```
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
#SBATCH --cpus-per-task=18
#SBATCH --ntasks=1


module load trinity	
module load trinity
#set in directory to where trimmed reads are stored
indir=/pine/scr/o/m/omtorano/exports/trimmed_reads
RUN=${SLURM_ARRAY_TASK_ID}
input=`ls ${indir}/*R1*gz | awk -F 'R1' '{print $1}'| sed -n ${RUN}p`
a=`basename ${input}`
#make unique out directory for each sample by amending sample name to trinity outdir !!out directory must contain 'trinity'!!
outdir=/pine/scr/o/m/omtorano/exports/assembly/${a}trinity/
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
```
The following trinity command is modeled after Natalie Cohen and Logan Whitehouse's examples. See the trinity github for details on specific flags. The grid_exec command will launch parallel jobs during the final phase of trinity and the .conf file has been optimized by Sara Haines. The locations of the .pl and .conf files (also included in this repo) should not change, so the paths to these files do not need to be modified unless they are moved. Within the trinity application there is the script TrinityStats.pl which will give details about the quality of the assembly (N50, # of contigs, longest contig, etc).
```
Trinity \
	--seqType fq \
	--max_memory 500G \
	--left ${input}R1*gz \
	--right ${input}R2*gz \
	--CPU ${SLURM_CPUS_PER_TASK} \
	--grid_exec "${griddir}/hpc_cmds_GridRunner.pl --grid_conf ${SLURM_SUBMIT_DIR}/hpc_conf/trinity_pipeline_longleaf_grid.conf -c" \
	--full_cleanup \
	--NO_SEQTK \
	--output ${outdir}
	/nas/longleaf/apps/trinity/2.8.6/trinityrnaseq-2.8.6/util/TrinityStats.pl ${outdir}/Trinity.fasta
```
Ideally Trinity will complete after ~2ish - 6ish days, even with the 'mail-user' sbatch option slurm will not send an email after job completion. Use the 'squeue -u <onyen>' command to check on the status of the array job. If there are hundreds of jobs in your queue this means trinity is running grid runner and is in phase 2. The .out files created by trinity also output the current status of the run. These files get huge as every single sequence is normalized and assembled. Errors during trinity runs can be cryptic, if a run runs out of time is will not say that in the .out file but it will in the .err file. If something is wrong with part of the command it will show up in the .out file but not the .err file. As written here if the job immidiatly fails it is likely an error in a path file. If it fails after running for days it could either be insufficent resources or a gridrunner error. For some reason the management of the gridrunner sometimes goes awry, you will be able to see if this has happened if there is a recursive_trinity.cmds.hpc-cache_success.__failures file in your trinity directory. As far as I can tell this is unavoidable and the huge difference in time to completion between trinity runs with and without the gridrunner flag makes using gridrunner worth it. Regardless of the source of the error, if the job has been running for days and fails the best option is to resubmit, triniity can pick up where it left off. This is handled by the next part of the script.
	
If the final trinity output file does not exist this loop will create a new script called trinity_rerun.bash that is the same as what was just run but with fewer resources (after 6 days most of trinity is complete). This script will then be submitted, after 'then' is 'sbatch ${a}trinity_rerun.bash'. Within the created bash script is an additional loop that again checks for the final trinity output, if it is not present it resubmits the same rerun script. The danger with this is that sometimes gridrunner errors can not be solved by resubmission. If the rerun script has run 3+ times cancel it, edit the rerun script by deleting the --grid_exec command and resubmit. In both the outside and rerun loops if the final trinity.fasta file does exist it will delete the intermediate files. This is important because trinity creates a huge number of intermediate files. 
```
 echo "Checking completion ..."
if [ ! -f ${outdir}.Trinity.fasta ]
then

cat<<EOF>${a}trinity_rerun.bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -J {a}trinity_
#SBATCH -o {a}trinity_.%j.out
#SBATCH -e {a}trinity_.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

Trinity \
	--seqType fq \
	--max_memory 100G \
	--left ${input}R1*gz \
	--right ${input}R2*gz \
	--CPU ${SLURM_CPUS_PER_TASK} \
	--grid_exec "${griddir}/hpc_cmds_GridRunner.pl --grid_conf ${SLURM_SUBMIT_DIR}/trinity_pipeline_longleaf_grid.conf -c" \
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
	/nas/longleaf/apps/trinity/2.8.6/trinityrnaseq-2.8.6/util/TrinityStats.pl ${outdir}/Trinity.fasta
	echo "I want to rm -r ${outdir}, but I'm scared"
fi


EOF

sbatch ${a}trinity_rerun.bash
else
    echo " ${outdir}.Trinity.fasta complete... Delete Trinity directory."
	/nas/longleaf/apps/trinity/2.8.6/trinityrnaseq-2.8.6/util/TrinityStats.pl ${outdir}/Trinity.fasta
	echo "I want to rm -r ${outdir}, but I'm scared"
fi
```
	
# CD hit
https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide
After assembling individual samples use CD hit to create "grand" assemblies, for Exports high yield created depth 1 and depth 4 assemblies. Along with combining indivisual assemblies into depth assemblies, this reduces contig redundancy.
```
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
```

# Annotate
Annotate with database of choice, in this project used KEGG for functional annotation and phylodb for taxonomic annotation. 
Make database of choice unto diamond database, use diamond to blast genes against database.

```
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
```
For this project then used keggannot and fastannotation.py to format annotation files.

# Align
Align individual samples with "grand" assemblies. Used salmon for this project, create "grand" assembly index, align with salmon quant. This script would have to be modified depending on the identifier within the sample name that dictates "group". For this project I am using 1 & 4 within the sample name to differentiate between depth 1 and 4 but this is project dependent. 
```
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
```

# Differential expression
This can be done numerous ways using different tools. For example, there are tools within trinity that run r scipts to get differential expression. See this tutorial for more information https://southgreenplatform.github.io/trainings/trinityTrinotate/TP-trinity/. Most differential expression is done in R using tools within the biocmanager package which is part of the bioconductor project https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html. At this point the route you take is so project dependent it is difficult to generalize. It is important to know what question you are trying to answer before blindly running scripts. For the rest of this tutorial I will only get into the beginning parts of using biocmanager. Here I use tximport to get all quant.sf files into correct format. To run r scripts on longleaf you must create a bash file that calls an r script. 
	
This simply gives the slurm options for the tximport.r file to be run with

```
#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-03:00:00
#SBATCH --mem=300G
#SBATCH --ntasks=3
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=
#SBATCH -J tximport
#SBATCH -o tximport.%A.out
#SBATCH -e tximport.%A.err

module load r/4.0.1
cd /pine/scr/o/m/omtorano/alignment/trinity/
outdir=/pine/scr/o/m/omtorano/tximport/

echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi



Rscript tximport.r
```
	
Within tximport.r are the actual commands being run. This loads the packages necessary for biocmanager and the data manipulation to follow.
```
BiocManager::install('tximport')
library(tximport)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(stringr)
	
```
Below the salmon output quant.sf files are being called and imported into r
```
samples<-list.files(path="//pine/scr/o/m/omtorano/alignment/trinity/", full.names=T)
files<-file.path(samples,"quant.sf")
names(files)<-str_replace(samples, "/pine/scr/o/m/omtorano/alignment/trinity/${input1}_quants_trinity","")%>%str_replace(".salmon","")
tx2gene <- read.delim("/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/kegg_annotations/mega_kegg.tsv")
txi.keggannot.transcriptlevel<-tximport(files,type="salmon",tx2gene=tx2gene[,c('query', 'KO')], txOut = TRUE)
```
The timing required for import varies by the number of quant.sf files and the size of these files, it is important to save intermediate steps, otherwise nothing will be saved until the final output and this script must be rerun fully to experiment with output.
```
txi.keggannot.transcriptlevel.1.rds <- saveRDS(txi.keggannot.transcriptlevel, "txi.keggannot.transcriptlevel.1.rds")
kegg_tl<-readRDS('txi.keggannot.transcriptlevel.1.rds')
rownames(kegg_tl) <- kegg_tl$X
colnames(kegg_tl)[colnames(kegg_tl) == 'X'] <-'query'
#phylodb annotations
tx2gene1<-read.delim("/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/phylodb_annotations/megaphylodb.tsv")
#import
txi.phylodb.transcriptlevel<-tximport(files, type="salmon", tx2gene =tx2gene1[,c("TrinityID", "Organism")],txOut=TRUE)
#save r data
txi.phylodb.transcriptlevel.1.rds<-saveRDS(txi.phylodb.transcriptlevel, "txi.phylodb.transcriptlevel.1.rds")
```
This was an experiment with a subset of the total phylodb dataset and can be modified and expreimented with
```								
phylo<-read.csv('sub.phylo.csv', header=T)
kegg<-read.csv('sub.phylo.csv', header=T)
rownames(kegg)<-kegg$X
colnames(kegg)[colnames(kegg)=='X']<-'TrinityID'

p<-read.delim('/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/phylodb_annotations/megaphylodb.tsv')
k<-read.delim('/proj/marchlab/projects/PUPCYCLE_2019/mega_annotations/annotations/kegg_annotations/mega_kegg.tsv')
colnames(keggannot)[colnames(keggannot)=='query']<-'TrinityID'
tpmphy<-merge(kegg, p, by ='TrinityID')
tpmphy<-merge(tpmphy,k, by ='TrinityID')
```
See this site http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data for more information about differential expression and biocmanager.
	

 More info on Spades and compiled useful unix commands to come :) 



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
  - ?

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
For trimming I used the tool trim_galore https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
This tool removes low quality reads and auto detects and removes common adapters. Here I will explain the trimgalore.bash script, to download a version that can be run (with minor edits to file path names etc.) see files listed above or git clone. Copy and pasting individual sections of code will not work correctly, this must be run as a script.

1) Longleaf is a big computer with lots of nodes, when you log on to longleaf you, and everyone else, are on the login node. Anything you do on the login node takes up computing resources, if longleaf is being really slow it could be because there are tons of people using it. Due to the limited resources, anything more than really resource light commands (i.e. moving files) should be done as a 'job'. A 'job' is basically what you call the commands you have saved in a file called a script once you ask it to be run. The way you ask for longleaf to run a job is by using the 'sbatch' command followed by the name of your script. Sbatch is a command that communicates with SLURM, the job scheduler that is used on longleaf. SLURM is just one of a few scheduling tools used by high performance computing platforms, but it is common enough that there are a lot of resources online. SLURM basically manages the logistics of allocating compute time, memory, nodes etc. See https://its.unc.edu/research-computing/techdocs/getting-started-on-longleaf/#Job%20Submission. To tell SLURM how many resources you need the following paragraph goes at the beginning of a script:
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
While '#' normally indicates lines within code that are not read, that isn’t the case here. The first line is called a shebang and for this script tells longleaf this is a bash script, other scripts (python, etc.) might look different. Below this are the SLURM options, here it's requesting a general node (instead of bigmem, snp, etc.), 1 node, 4 hours of time, 10 gigs of memory, 5 tasks, indicating that we are requesting an array and giving the length of the array (here I have 21 samples), assigning the job name to be trim, assigning the name of the '.out' file to be trim.jobnumber_arraynumber.out, same for the error file but .err, requesting longleaf send an email when the job starts, ends, or fails, and give the email address to send the notification to. Most (all?) functions when run create a 'output stream' and 'error stream', basically any text output from a function working properly will be part of the 'output stream', and errors will be part of the 'error stream'. Instead of printing to the terminal screen these streams are passed into files that will appear in your working directory. This job creates 21 .out and 21 .err files, hence it is important to add the unique job number to the name. 

2) load necessary modules for trim_galore, auto loads python & cutadapt
```
module load trim_galore
module load pigz
```
Longleaf has many software tools available that can be used by adding them to your working space with 'module load' (module add does the same thing). Other useful module management commands: module avail - lists allllll available modules (modify by doing module avail r* etc. to search for modules starting with r etc.), module list - shows all the modules you have loaded in your working space (loading a module in a job does not mean it will be loaded in your login session), module rm <module name> - removes specific module, module purge - removes all modules loaded in your space. Most longleaf modules are stored in longleaf here: /nas/longleaf/apps/

3) set 'in' directory to where raw reads are, if you need to transfer reads from /proj run rsync -r (recursive, needed if transferring all files in directory) /from/path /to/path
```
rsync -r /proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/Reads /pine/scr/o/m/omtorano
indir=/pine/scr/o/m/omtorano/exports/reads
outdir=/pine/scr/o/m/omtorano/exports/trimmed_reads
```
4) This creates the out directory if it does not already exist. Echo is a Linux command which basically means print - i.e. this line will print in your .out file "checking if out directory exists", this is not necessary for creating the directory, just helpful to tell you what it’s doing. Below this is a for loop and conditional saying if out directory does not exist make it, if it does exist print "... exists". The"!" basically means 'not this', in this case if outdir does not exist. 
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
5) Set variable called 'RUN to = slurm array task id (https://slurm.schedmd.com/job_array.html) which is a very cool tool for submitting a bunch of jobs at once, more explanation below. This isn’t really necessary but makes the following line a bit cleaner.
```
RUN=${SLURM_ARRAY_TASK_ID}
```
6) Parse and defines input as the sample name that will be input to trim_galore. ls ${indir}/&ast;R1&ast; will list all files in in directory that have R1 anywhere in the name, the * wildcard is used here to mean find R1 embedded in any string of characters before or after. &ast;R1 would only look for R1 occurring at the end of any string of characters, R1&ast; would only look for R1 at the beginning. This is done to get each sample name listed one time (there will be an R1 and R2 in the indir assuming paired end reads). The '|' character 'pipes' commands, the output of ls is given to awk. awk here is being used to cut path file name at R1, -F is the field separator which is set here to R1. This results in names being split before and after the R1, '{print $1}' prints the first element of the split name, now we have each unique sample ID listed one time. sed here is being used to isolate the 'changed' input names. The slurm array task id will basically run through a list of all of the unique sample names, to input them one at a time sed -n followed by p prints only the file being processed by the slurm array. sed and awk are useful Linux tools however there are multiple ways to accomplish this string parsing to get the unique sample names.
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
If working on one project in scratch it would not be necessary to make further subdirectories for reads as is written here (reads is a subdirectory of exports in this example).

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
set out directory to where fastqc output should go, can make 2 directories to keep raw and trimmed output separate
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
run fastqc on untrimmed raw reads and trimmed reads, reads will run simultaneously, set threads to # of reads
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
The multiqc step will result in an .html file that can be dragged and dropped into a browser. This file contains the compiled fastq reports from all trimmed and raw reads. If there are certain samples that have low quality or vary from the rest of the samples this may indicate that these sample need to be removed from further analysis. Full multiqc output for the High Yield 2020 samples is saved in HighYield2020_multiqc_report.html, below is an example of what it looks like.. 
	![tempsnip](https://user-images.githubusercontent.com/48129653/123166467-71ad1d00-d443-11eb-9bf7-d6b0f8355739.png)


# Assembly with Trinity
The goal of assembly in this work flow is to create a de novo assembly from all samples in a 'group' which you can map back to in the alignment step. For the Exports High Yield 2020 samples were collected on different days from two different depths, 1 & 4. I am interested in comparing between days, not depths, I am therefore creating an assembly of all depth 1 samples and all depth 4 samples. This way I can compare samples collected at a single depth, instead of between both depths. 

Trinity is one option for assembly tools, the other I will go through here (next section) is Spades. Trinity is a multipurpose tool that can do a lot more than assembly, read more here https://github.com/trinityrnaseq/trinityrnaseq/wiki. During assembly it runs through different phases (see wiki) that can be split into different jobs depending on how big an assembly you are trying to make. For this project the stats of my assemblies are below:

![image](https://user-images.githubusercontent.com/48129653/119878716-7809b100-bef8-11eb-97c4-d12eb5d70787.png)

Totals were calculated based on averages (total sequence length = average sequence length * number of samples), averages were calculated during multiqc step.

	



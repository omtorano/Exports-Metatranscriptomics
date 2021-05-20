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
		- You'll need the folloing info, -> are my suggestions
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
	- There are many ways to do this. Some options include using longleaf desktop https://ondemand.rc.unc.edu/pun/sys/dashboard/ (Sarah's preferred method), doing everything (navigating directories, editing files, etc.) from a terminal (Johnson's preffered method), doing hybrid terminal and file manager (my preferred method, on Windows I use a GitBash and WinSCP file manager, it works for me). 
- Make a new directory within /proj/marchlab/projects for your project, or new subdirectory within one of the existing project directories. Then download sequences if sequencing was done by an outside service (not UNC's High Throughput Sequencing facility HTSF)
	- For example - when I got sequences for the High Yield 2020 EXPORTS project I made a directory within the existing Exports project folder called /HighYield2020. The sequencing was completed at Genewiz, who sent email detailing how to transfer files with sftp. This is an example of what this would looked like (JUST AN EXAMPLE, DO NOT RUN):
```
### log into longleaf

ssh omtorano@longleaf.unc.edu

### password prompt will come up automatically, you will not be able to see characters as you type them

Password: <enter password here>

### cd changes directory

cd /proj/marchlab/projects/EXPORTS/metatranscriptomics/

### mkdir makes a directory, here I am making the directory HighYield2020 and another directory within that called Reads on the following line

mkdir HighYield2020
mkdir HighYield2020/Reads

### the beginning of the directions from Genewiz - connect to their server via sftp

sftp a_marchetti_lab_gmail@sftp.genewiz.com

### password prompt will come up automatically, enter the password they sent via email

password: <enter password genewiz sent here>

### lcd sets the local directory (where you want your reads to go), does not change working directory (does not change where you are)

lcd /proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/Reads

### ls lists items in the current working directory, in this case I just logged onto Genewiz so I am seeing the list of directories
### I have access to - likely there will be one and it will have the name of the Genewiz project number

ls

### change directory into what ever directory was just listed by ls

cd <genewiz project folder name that just came up from ls command>

### mget transfers files from working directory to whatever you set as your local directory. The '*' is a Linux wildcard that means 'get all the stuff'

mget * #I dont remember if there is a space here

### end the sftp connection to transfer yourself back to longleaf

quit #Is this the command that ends sftp connection?
```
## Getting Started
Rules for the Marchetti Lab /proj space
- Do not work in the proj space, it is only for storage of 'final product' files from each stage in pipeline. We have a history of running out of space, storage availability can be tracked here https://rc-storage-info.its.unc.edu:32000/home.
- Instead work in scratch (/pine/scr/<o>/<n>/<onyen>) or home directory (/nas/longleaf/home/<onyen>). Scratch has tons of space but inactive files will be deleted by ITS after ~30days, they send an email before deletion.
- Use the rsync command to move files, mv deletes files from their original location and if interupted will result in data loss.
- Use the rm (remove) command CAREFULLY, permanantly deletes files.
- The recommended directory organization is to have directory for the project and seperate sub directories for code, reads, assemblies, annotations, and alignment (optional trimmed reads and fastqc output).
- See final section of page for useful stuff.

# Trimming
For trimming I used the tool trim_galore https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
This tool removes low quality reads and auto detects and removes common adapters. Here I will explain the trimgalore.bash script, to download a version that can be run (with minor edits to file path names etc.) see files listed above or git clone. 

Longleaf is a big computer with lots of nodes, when you log on to longleaf you, and everyone else, are on the login node. Anything you do on the login node takes up computing resources, if longleaf is being really slow it could be because there are tons of people using it. Due to the limited resources, anything more than really resource light commands (ie. moving files) should be done as a 'job'. A 'job' is basically what you call the commands you have saved in a file called a script once you ask it to be run. The way you ask for longleaf to run a job is by using the 'sbatch' command followed by the name of your script. Sbatch is a command that communicates with SLURM, the job scheduler that is used on longleaf. SLURM is just one of a few scheduling tools used by high performance computing platforms, but it is common enough that there are a lot of resources online. SLURM basically manages the logistics of allocating compute time, memory, nodes etc. See https://its.unc.edu/research-computing/techdocs/getting-started-on-longleaf/#Job%20Submission. To tell SLURM how many resources you need the following paragraph goes at the beginning of a script:


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

While '#' normally indicates lines within code that are not read, that isnt the case here. The first line is called a shebang and for this script tells longleaf this is a bash script, other scripts (python, etc.) might look different. Below this are the SLURM options, here it's requesting a general node (instead of bigmem, snp, etc.), 1 node, 4 hours of time, 10 gigs of memory, 5 tasks, indicating that we are requesting an array and giving the length of the array (here I have 21 samples), assigning the job name to be trim, assigning the name of the '.out' file to be trim.jobnumber_arraynumber.out, same for the error file but .err, requesting longleaf send an email when the job starts, ends, or fails, and give the email address to send the notification to. Most (all?) functions when run create a 'output stream' and 'error stream', basically any text output from a function working properly will be part of the 'output stream', and errors will be part of the 'error stream'. Instead of printing to the terminal screen


### load necessary modules for trim_galore, auto loads python, cutadapt

module load trim_galore
module load pigz


#### set in directory to where raw reads are, if you need to transfer reads from /proj run line 94, rsync -r (recursive, needed
### if transfering all files in directory) /from/path /to/path

#rsync -r /proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/Reads /pine/scr/o/m/omtorano

indir=/pine/scr/o/m/omtorano/exports/reads
outdir=/pine/scr/o/m/omtorano/exports/trimmed_reads

### lines 129-1135 create out directory if it does not already exist, echo is a linux command which basically means print - ie 
### line 104 will print in your .out file "checking if out directory exists". 105-111 is a for loop and conditional saying if 
### out directory does not exist make it, if it does exist print "... exists". 

echo "Checking if ${outdir} exists ..."
if [ ! -d ${outdir} ]
then
    echo "Create directory ... ${outdir}"
    mkdir -p ${outdir}
else
    echo " ... exists"
fi

### setting variable called 'RUN to = slurm array task id (info here https://slurm.schedmd.com/job_array.html) which is a very cool tool for submitting a bunch of jobs at once

RUN=${SLURM_ARRAY_TASK_ID}

### ls ${indir}/*R1* will list all files in in directory that have R1 anywhere in the name, the * wildcard is used here to mean 'find R1 embedded in any string of characters 
### before or after'. *R1 would only look for R1 occuring at the end of any string of characters, R1* would only look for R1 at the beginning. The '|' character 'pipes' 
### commands, the output of ls is given to awk. awk here is being used to cut path file name at R1, -F is the field separator which is set here to R1, which results in names
### being split before and after the R1, '{print $1}' prints the first element of the split name. sed here is being used to isolate the 'changed' input names. The slurm array 
### task id will basically run through a list of all of the files in the in directory, to input them one at a time sed -n followed by p prints only the file being processed by 
### the slurm array. TBH I dont completely understand the nuance but it works. 

input=`ls ${indir}/*R1* | awk -F 'R1' '{print $1}'| sed -n ${RUN}p`

### The -e flag here allows echo to interprate the backslash escapes, here used to print on a new line "\n", other uses include "\t" if you want to echo something but have it
### be tab seperated

echo -e "\nRun ID: ${RUN}"
echo -e "\nSample: ${input}"

### the actual trimming part, -j tells it how many cores to use, I think due to python versions it will actually only run one core as written here

trim_galore -j 4 \
	--paired ${input}R1* ${input}R2* \
	-o  ${outdir}
	

```
text  

# dkfjlj  
## header!!!  
###### so much header  

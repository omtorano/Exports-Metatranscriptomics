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
#log into longleaf
ssh omtorano@longleaf.unc.edu
#password prompt will come up automatically, you will not be able to see characters as you type them
Password: <enter password here>
#cd changes directory
cd /proj/marchlab/projects/EXPORTS/metatranscriptomics/
#mkdir makes a directory, here I am making the directory HighYield2020 and another directory within that called Reads on the following line
mkdir HighYield2020
mkdir HighYield2020/Reads
#the beginning of the directions from Genewiz - connect to their server via sftp
sftp a_marchetti_lab_gmail@sftp.genewiz.com
#password prompt will come up automatically, enter the password they sent via email
password: <enter password genewiz sent here>
#lcd sets the local directory aka where you want your reads to go, does not change working directory (aka does not change where you are)
lcd /proj/marchlab/projects/EXPORTS/metatranscriptomics/HighYield2020/Reads
#ls lists items in the current working directory, in this case I just logged onto Genewiz so I am seeing the list of directories I have access to - likely there will be one and it will have the name of the Genewiz project number
ls
#change directory into what ever directory was just listed by ls
cd <genewiz project folder name that just came up from ls command>
#mget transfers files from working directory to whatever you set as your local directory. The '*' is a Linux wildcard that means 'get all the stuff'
mget * #I dont remember if there is a space here
#end the sftp connection to transfer yourself back to longleaf
quit #Is this the command that ends sftp connection?
```


# Trimming
For trimming I used the tool trim_galore https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
This tool removes low quality reads and auto detects and removes common adapters.

```
#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=00-4:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=5
#SBATCH --array=1-21 #array #s needed?
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

#create ourdirectory
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
input=`ls ${indir}/*R1* | awk -F 'R1' '{print $1}'| sed -n ${RUN}p`


echo -e "\nRun ID: ${RUN}"
echo -e "\nSample: ${input}"


trim_galore -j 4 \
	--paired ${input}R1* ${input}R2* \
	-o  ${outdir}
	

```
text  

# dkfjlj  
## header!!!  
###### so much header  

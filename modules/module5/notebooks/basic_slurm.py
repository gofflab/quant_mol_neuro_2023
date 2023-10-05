#%% [markdown]
# # Module 5 - Read QC with a helpful hand from slurm
# This notebook will introduce you to the basics of using slurm to submit jobs to the cluster.  We will use slurm to run the read qc program `fastqc` on a collection of fastq files provided for you.
#
# ## Learning Objectives
# - Familiarize yourself with the basic slurm commands
# - Learn how to submit jobs to the cluster using the slurm command `sbatch`
# - Learn how to monitor jobs using the slurm command `squeue`
# - Understand where and how to find the output of your jobs.
# - Gain experience in reviewing the output of fastqc to assess the quality of the reads in a fastq file.
#

#%% [markdown]
# ## The Task
# We have provided you with a set of fastq files in the directory `/data/me440_lgoff2/datasets/` on the rockfish cluster.  Your task is to copy the fastq files to your working directory and run `fastqc` on each of these files and review the output to assess the quality of the reads in each file.  
# 
# You will submit these jobs to the cluster using slurm.
#
# ## The Tools
# We will be using the following tools:
# - `sbatch` - submit a job to the cluster (slurm is already provided/installed on the cluster so we do not need to add anything to our environment)
# - `fastqc` - a program that assesses the quality of reads in a fastq file
# - `multiqc` - a program that aggregates the output of multiple steps in a bioinformatics workflow (including fastqc runs) into a single report.
#
# To install `fastqc` and `multiqc` on your rockfish account, we will use the following command:
#
# _You will only need to run this once to install._

#%%
%%bash
mamba install -c bioconda fastqc multiqc

#%% [markdown]
# _You will only need to run this once to install._

#%% [markdown]
# ## The Data
# Next, Let's make a copy of the fastq files which are already stored on the cluster in a shared directory.
#
# We'll start by creating a new directory 'data' in our current working directory to store the fastq files.
# 

#%%
%%bash
mkdir data

#%% [markdown]
# Now, we'll copy the fastq files from the shared directory to our local directory.
#
# These are the raw RNA-Seq reads for the HippoSeq dataset [GSE74985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74985).
# 
# The directory that contains the .fastq files is `/data/me440_lgoff2/datasets/RNA-Seq/data/raw/GSE74985`.  We can use the `cp` command to copy any file with a .fastq.gz extention from this directory to our new local directory.

#%%
%%bash
cp /data/me440_lgoff2/datasets/RNA-Seq/data/raw/GSE74985/*.fastq.gz data/

#%% [markdown]
# Let's take a look at the files we just copied over.

#%% 
%%bash
ls data/

#%% [markdown]
# For this exercise we want to run the read quality control program `fastqc` on each of these files.  
# 
# `fastqc` is a program that assesses the quality of reads in a fastq file.  It generates a report that can be used to explore the read quality.  
# 
# First let's create a directory to store the output of `fastqc` for each of the files.  We'll call this directory `fastqc_output`.

#%%
%%bash
mkdir fastqc_output

#%% [markdown]
# To generate the fastqc reports we need to run the following command for each file:
# 
# `fastqc -o fastqc_output <fastq_file>`
#
# Let's try the first file `SRR2916027.fastq.gz`

#%%
%%bash
fastqc -o fastqc_output data/SRR2916027.fastq.gz

#%% [markdown]
# This generates a report for the file `SRR2916027.fastq.gz` and stores it in the directory `fastqc_output` and takes ~3 minutes to run.
# 
# Let's take a look at the .html report that was generated. To do so, we can either download the file and open in a browser, or preview the file in vscode.  Let's try the latter.
#
# (Cmd+Shift+P) -> HTML: Open Preview

#%% [markdown]
# We still have a number of files for which to run this report.  We could run each job one at a time, but that would take a lot of babysitting.
#
# We could _also_ write a for loop to run each job in succession but it would still run each job sequentially.
#
# Instead, we can use slurm to submit each job to the cluster and let slurm manage the resources, scheduling, and execution of each job in parallel.

#%% [markdown]
## SLURM: An Introduction
#
# Slurm (Simple Linux Utility for Resource Management) is an open-source job scheduler that allocates compute resources on clusters for queued user jobs. 
# Slurm has become a standard for supercomputing environments, providing both resource management and job scheduling. Slurm is used on the [rockfish cluster](http://edulogin.arch.jhu.edu).
# 
# Slurm is a command-line tool that can be used to submit, monitor, and cancel jobs. 

# Slurm is _primarily_ useful when we need to run a large number of the same type of jobs that can be run independently of each other. This type of problem is called 'embarrasingly parallel'. Instead of running each of the jobs in succession, we can submit them all at once and let Slurm manage the resources, scheduling, and execution of each job in parallel.
# 
# This is the case for many bioinformatics pipelines, including for example, the alignment of many samples to a reference genome. 

#%% [markdown]
### 1. Submitting Jobs with `sbatch`
#
# `sbatch` is used to submit a job for later execution. 
# A 'job' is a set of commands that you would like to execute. 
#
# Typically, you provide `sbatch` with a script that provides directives and commands, but for simpler use cases, the `--wrap` argument allows direct command-line submission.
#
# **Submission using `--wrap`**:
#
# let's say we wanted to run the following command:
#
# `echo 'Hello World!'`
#
# We could submit this job using `sbatch` as follows:
#%% 
%%bash
sbatch --wrap="echo 'Hello World!'"

#%% [markdown]
# Notice the only information that we get from this is an acknowledgement that the job was submitted and a job id number.
#
# By default, the output of the job is written to a file called `slurm-<job_id>.out` in the current working directory. 
#
# Let's take a look at the contents of this 'output' file.

#%% [markdown]
# ### Back to the matter at hand
# Let's revisit our fastqc example.  We want to run the following command for each file:
# 
# `fastqc -o fastqc_output <fastq_file>`
#
# We can submit this job using `sbatch` as follows:
#
# `sbatch --wrap="fastqc -o fastqc_output <fastq_file>"`
#
# Let's rerun the first file `SRR2916027.fastq.gz` using `sbatch`:

#%%
%%bash
sbatch --wrap="fastqc -o fastqc_output data/SRR2916027.fastq.gz"

#%% [markdown]
# We can check the status of our jobs using `sacct`.  Let's take a look.

#%%
%%bash
sacct

#%% [markdown]
# We can see that our job is currently running.  We can also see the job id number, the user who submitted the job, the start time, the partition, the state, the exit code, and the elapsed time.

#%% [markdown]
# Ok, let's save ourselves some time and setup a `for` loop to submit each of the fastqc jobs for each of the `.fastq.gz` files in our `data` directory.
#

#%%
%%bash
for file in data/*.fastq.gz
do
	sbatch --wrap="fastqc -o fastqc_output $file"
done

#%% [markdown]
# Let's check the status of our jobs again using `sacct`.

#%%
%%bash
sacct

#%% [markdown]
### 2. Monitoring Jobs with `squeue`
#
#`squeue` displays the status of jobs. 
#
#**Basic usage**:
#
#%% [bash]
%%bash
squeue

#%%
%%bash
squeue


#%% [markdown]
# **Filter by User**:

#%% 
%%bash
squeue -u lgoff2


#%% [markdown]
# **Customize Display**:

#%% 
%%bash
squeue -o "%.4i %.9P %.30j %.8u %.2t %.10M %.6D %R"


#%%
### 3. Cancelling Jobs with `scancel`

`scancel` is used to terminate jobs.

**Cancel a Specific Job**:

#%% 
%%bash
scancel <job_id>


#%% [markdown]
# **Cancel All Jobs for a User**:

#%% 
%%bash
scancel -u lgoff2


#%%
### 4. Inspecting Cluster with `sinfo`
# `sinfo` provides an overview of Slurm nodes and partitions.
# **Basic Usage**:

#%% 
%%bash
sinfo


#%% [markdown]
# **Detailed View**:
#

#%% 
%%bash
sinfo -l


#%% [markdown]
# **Node-specific Information**:
#

#%% 
%%bash
sinfo -N -l


#%% [markdown]
### 5. Running Jobs Interactively with `srun`
#
# `srun` can be used to initiate job steps in real time.
#
# **Basic Usage**:

#%% 
%%bash
srun --pty /bin/bash


#%% [markdown]
# **Specify Resources**:

#%% 
%%bash
srun --partition=general --nodes=2 --ntasks-per-node=2 /bin/bash


#%% [markdown]
# ## MultiQC
# MultiQC is an exceptionally useful tool for aggregating the output of multiple steps or samples in a bioinformatics workflow (including fastqc runs) into a single report.
#
# MultiQC traverses a directory and searches for output from common bioinformatics tools, extracts the data, and generates a single report summarizing the data.
#
# Let's use MultiQC to aggregate the output of our fastqc runs into a single report.

#%%
%%bash
multiqc .

#%% [markdown]
# Let's take a look at the summary report that was generated.

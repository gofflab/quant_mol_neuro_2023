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

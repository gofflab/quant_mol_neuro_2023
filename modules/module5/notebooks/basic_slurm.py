
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

#%% 
%%bash
sbatch --job-name=wrapped_job --wrap="<bash command>"


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
squeue -u your_username


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
scancel job_id


#%% [markdown]
# **Cancel All Jobs for a User**:

#%% 
%%bash
scancel -u your_username


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

#%% [markdown]
# # Module 5 - Read QC with a helpful hand from slurm
#
# In this notebook, you will learn the basics of using SLURM (Simple Linux Utility for Resource Management) to submit jobs to the Rockfish cluster.
#
# **This notebook must be run on the Rockfish cluster (edulogin.arch.jhu.edu) using the login credentials provided on Monday.**
#
# **Running this notebook on your local machine is not be possible.**
#
# We will use SLURM to run the read quality control program `fastqc` on a collection of fastq files that we provide.
#
# ## Learning Objectives
# - Familiarize yourself with the basic SLURM commands.
# - Learn how to submit jobs to the cluster using the SLURM command `sbatch`.
# - Learn how to monitor jobs using the SLURM command `sacct`.
# - Understand where and how to find the output of your jobs.
# - Gain experience in reviewing the output of fastqc to assess the quality of the reads in a fastq file.

#%% [markdown]
# ## The Task
# We have provided you with a set of raw read fastq files in the directory `/data/me440_lgoff2/datasets/` on the Rockfish cluster.
# Your task is to copy the fastq files to your working directory, run `fastqc` on each of these files, and review the output to assess the quality of the reads in each file.
#
# You will submit these jobs to the cluster using SLURM.
#
# ## The Tools
# We will be using the following tools:
# - `sbatch` - submit a job to the cluster (SLURM is already provided/installed on the cluster so we do not need to add anything to our environment)
# - `fastqc` - a program that assesses the quality of reads in a fastq file
# - `multiqc` - a program that aggregates the output of multiple steps in a bioinformatics workflow (including fastqc runs) into a single report.
#
# To install `fastqc` and `multiqc` on your Rockfish account, we will use the following command:
#
# _You will only need to run this once to install._

#%%
%%bash
mamba install -c bioconda fastqc multiqc

#%% [markdown]
# **You will only need to run this once to install.**

#%% [markdown]
# ## The Data
# Next, we will copy the fastq files which are already stored on the cluster in a shared directory.
#
# Let's start by creating a new directory 'data' in our current working directory to store the fastq files.

#%%
%%bash
mkdir data

#%% [markdown]
# Now, let's copy the fastq files from the shared directory to our local directory.
#
# These are the raw RNA-Seq reads for the HippoSeq dataset [GSE74985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74985).
#
# The samples for this study were each sequenced once (1 run per sample) on the Illumina sequencing platform to generate single-end reads of 100bp in length.
#
# The directory that contains the gzip-compressed .fastq.gz files is `/data/me440_lgoff2/datasets/RNA-Seq/data/raw/GSE74985`.
#
# We'll use the `cp` command to copy all of the .fastq.gz files from the shared directory to our new local directory.

#%%
%%bash
cp /data/me440_lgoff2/datasets/RNA-Seq/data/raw/GSE74985/*.fastq.gz data/

#%% [markdown]
# Let's take a look at the files we just copied over.

#%% [bash]
%%bash
ls data/

#%%
%%bash
ls data/ | wc -l

#%% [markdown]
# There should be 24 *.fastq.gz files in the `data` directory.

#%% [markdown]
# For this exercise, we want to run the read quality control program `fastqc` on each of these files.
#
# `fastqc` is a program that parses the reads in a fastq file and outputs some summary statistics,
# metrics, and heuristics to tell us more about what is in the file, and the quality of the reads. It generates a report that can be used to identify systematic issues/errors that can be common in the library preparation or sequencing of the reads.
#
# Let's run fastqc to see if we can get an idea of the quality of the reads in one of the files.
#
# First, let's create a directory to store the output of `fastqc` for each of the files. We'll call this directory `fastqc_output

#%%
%%bash
mkdir fastqc_output

#%% [markdown]
# To generate the fastqc reports we need to run the following command for each file:
#
# ```
# fastqc -o fastqc_output <fastq_file>
# ```
#
# The `-o` flag tells `fastqc` where to store the output files, in this case,the directory `fastqc_output`.
#
# _Remember, for most tools we can either use the `man` command or the `--help` flag to get more information about the tool and how to use it.
# `fastqc` is no exception here and if we want to find other arguments that are available we can run:
#
# ```
# fastqc --help
# ```
#
# Let's try the first file `SRR2916027.fastq.gz`

#%%
%%bash
fastqc -o fastqc_output data/SRR2916027.fastq.gz

#%% [markdown]
# This generates a report for the file `SRR2916027.fastq.gz` and stores it in the directory `fastqc_output`. This process takes approximately 3 minutes to run.
#
# Let's take a look at the .html report that was generated. To do so, we can either download the file and open in a browser, or preview the file in vscode. Let's try the latter.
#
# (Cmd+Shift+P) -> HTML: Open Preview (_You may have to install the `HTML Preview` extension for VSCode_)

#%% [markdown]
# We still have a number of files for which to run this report.
# We could run each job one at a time, which would require constant monitoring.
#
# Alternatively, we could write a for loop to run each job in succession, but it would still run each job sequentially (~3min * 24 jobs = ~72 minutes).
#
# Instead, we can use SLURM to submit each job to the cluster and let SLURM manage the resources, scheduling, and execution of each job in parallel.

#%% [markdown]
# ### SLURM: An Introduction
#
# [Slurm](https://slurm.schedmd.com/quickstart.html) (Simple Linux Utility for Resource Management) is an open-source job scheduler that allocates compute resources on clusters for queued user jobs.
# Slurm has become a standard for supercomputing environments, providing both resource management and job scheduling. Slurm is used on the [rockfish cluster](http://edulogin.arch.jhu.edu).
#
# Slurm is a command-line tool that can be used to submit, monitor, and cancel jobs.
#
# Slurm is _primarily_ useful when we need to run a large number of the same type of jobs that can be run independently of each other.
# This type of problem is called 'embarrasingly parallel'. Instead of running each of the jobs in succession, we can submit them all at once and let Slurm manage the resources, scheduling, and execution of each job in parallel.
#
# This is the case for many bioinformatics pipelines, including for example, the alignment of many samples to a reference genome.

#%% [markdown]
# #### Job Submission with `sbatch`
#
# The `sbatch` command is used to submit a job for execution to the SLURM scheduler.
# A 'job' in this context is a set of commands you wish to execute.
#
# You typically provide `sbatch` with a script containing directives and commands.
# However, for simpler cases like ours, the `--wrap` argument allows direct submission of command-line calls.

# #### Submission using `--wrap`:
#
# Let's say we want to run the following command:
#
# `echo 'Hello World!'`
#
# We can submit this job using `sbatch` as follows:
#%%
%%bash
sbatch --wrap="echo 'Hello World!'"

#%% [markdown]
# Upon submission, SLURM returns a `job id` number. The output of the job is written to a file named `slurm-<job_id>.out` in the current working directory. This file captures the STDOUT, which is what would normally be printed to the terminal. In our case, it contains the string 'Hello World!'.
#
# Let's examine the contents of this 'output' file.

#%% [markdown]
# ### Applying `sbatch` to the `fastqc` example
# We want to run the following command for each file:
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
# The job is now 'submitted' to the SLURM scheduler and we have a new `job_id` number.
# We can check the corresponding `.out` file to see the output (STDOUT) that would have been printed to the screen.
#
# Note that this job execution doesn't tie up our terminal resources as it runs elsewhere on the cluster,
# allowing us to continue working.

#%% [markdown]
# To check the status of our jobs
# we use the SLURM command `sacct` ([SLURM 'Accounting'](https://slurm.schedmd.com/sacct.html)).

#%%
%%bash
sacct

#%% [markdown]
# We can see that our job is currently running.
# The output includes the job id number, the user who submitted the job, the start time, the partition,
# the state, and any exit codes (errors) that the job may have produced (normal exit is 0).

#%% [markdown]
# To streamline the submission and execution of the `fastqc` jobs for the `.fastq.gz` files in our `data` directory,
# we can use a bash `for` loop to iterate over each file and submit a job.

#%%
%%bash
for file in data/*.fastq.gz
do
	sbatch --wrap="fastqc -o fastqc_output $file"
done

#%%
# Python version

import subprocess
from pathlib import Path

for file in Path("data").iterdir():
    if file.suffix == ".fastq.gz":
        subprocess.run(
            f'sbatch --wrap="fastqc -o fastqc_output {file}"',
            check=True,
            shell=True
		)

#%% [markdown]
# Let's check the status of our jobs again using `sacct`.
#%%
%%bash
sacct

#%% [markdown]
# Now we can see more jobs in the queue.
# Some are running, some might have completed, and some might still be 'PENDING'.
# SLURM balances the available resources on the cluster to ensure all jobs complete in a timely manner,
# queuing some jobs until resources become available.

# SLURM generally follows a 'first-come, first-served' model,
# where jobs submitted first have higher priority for execution.

#%% [markdown]
# ### Monitoring All Jobs with `squeue`
#
# While `sacct` is useful for monitoring your own jobs,
# `squeue` displays the status of all 'active' jobs on the cluster, including those submitted by other users.

#%% [bash]
%%bash
squeue

#%% [markdown]
# The output shows a list of jobs in the queue,
# distributed across different compute nodes (NODELIST) on the cluster.
# This gives an idea of the cluster load and available resources.
# To focus on _your_ specific jobs, use the `-u` flag to filter by user.

#%%
%%bash
squeue -u lgoff2

#%% [markdown]
#### Cancelling Jobs with `scancel`
#
#`scancel` is used to terminate 'RUNNING' or 'PENDING' jobs.
#
# The `scancel` command is used to terminate 'RUNNING' or 'PENDING' jobs.

##### Cancel a Specific Job:
# To cancel a specific job, provide the job id number to the `scancel` command:

# ```
# scancel <job_id>
# ```

#%% [markdown]
##### Cancel All Jobs for a User:
# If you've made an error while submitting a large number of jobs, you can cancel all jobs for a specific user:
#
# ```
# scancel -u lgoff2
# ```
#
# However, caution is advised as this will cancel **all** of your jobs, including those that are maintaining your ssh tunnel to allow VSCode to connect. Doing this will disconnect you from the cluster and you will have to reconnect.

#%% [markdown]
#### Inspecting Cluster with `sinfo`
# As we delve deeper into using the cluster and SLURM,
# it may be useful to understand more about the cluster's configuration and the resources available for SLURM.
#
# The `sinfo` command provides an overview of SLURM nodes (compute nodes) and partitions ('queues' for job submission).

#%%
%%bash
sinfo

#%% [markdown]
# The output tells us that there is one (unique) partition or 'queue' for job submission called `defq` (for 'default queue').
# It is currently available (up) with a 2-hour time limit per job (`2:00:00`).
#
# It also informs us about the number of compute nodes available for this partition (6)
# and the 'NODELIST' provides the 'names' of the compute nodes included in this partition.


#%% [markdown]
###### Node-specific Information
# We can also get more information about the compute nodes using the `-N` flag
#
# `sinfo -N`
#
# Let's also add the `-l` flag to give us more information in a 'long' format.

#%%
%%bash
sinfo -N -l

#%% [markdown]
# Here we can see that each of the compute nodes has 24 available CPUs and 91552 MB (~92 Gb) of memory available for jobs.
#
# These are the resources that slurm is managing and allocating to jobs that are submitted to the scheduler.

#%% [markdown]
# ### Running Jobs Interactively with SLURM

# We've already seen `srun` in action when we used it to create an interactive session and the SSH tunnel on the cluster. Let's revisit this command to understand its components better.

# ```bash
# ssh lgoff2@edulogin.arch.jhu.edu "srun --time=2:00:00 --mem-per-cpu=4GB --cpus-per-task=2 VSCode-linux-x64/bin/code tunnel --accept-server-license-terms"
# ````
#
# The above command was used to ssh into the cluster using my credentials (lgoff2@edulogin.arch.jhu.edu), and then immediately execute a call to srun on the login node.

# Here, srun initiates an interactive session on the cluster with the following resources:

# - --time=2:00:00 - 2 hours of walltime
# - --mem-per-cpu=4GB - 4 GB of memory per CPU
# - --cpus-per-task=2 - 2 CPUs per task
#
# By default, `srun`` starts an interactive session and runs the provided commands.
# In this instance, we're asking srun to create an interactive session and run the following command immediately:
#
# ```
# VSCode-linux-x64/bin/code tunnel --accept-server-license-terms
# ```

#%% [markdown]
# This command starts the VSCode server on the cluster and creates the SSH tunnel, enabling us to connect to the server.
#
# Interactive sessions on the cluster can be helpful when debugging or developing job submissions.
#

#%% [markdown]
# ## MultiQC
# Having successfully executed fastq on each of our fastq files and generated a report for each file,
# we'd like to examine these reports to identify any unusual patterns or suspect read quality in these samples.

# However, parsing through 24 .html files to find patterns/trends can be tedious.
# This is where MultiQC comes in handy. [MultiQC](https://multiqc.info/) is a powerful tool for aggregating the output of multiple steps or samples in a bioinformatics workflow (including fastqc runs) into a single report.

# MultiQC traverses a directory, searches for output from common bioinformatics tools, extracts the data, and generates a single report summarizing the data.
# Let's use MultiQC to aggregate the output of our fastqc runs into a single report, which will simplify comparisons across our samples.
#
#%%
%%bash
multiqc .

#%% [markdown]
# Let's take a look at the summary report that was generated.
#
# See anything useful or interesting about the samples?

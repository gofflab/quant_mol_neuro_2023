#%% [markdown]
# # Module 5 - Read QC with a helpful hand from slurm
# 
# This notebook will introduce you to the basics of using slurm to submit jobs to the cluster. 
# 
# __This notebook must be run on the Rockfish (edulogin.arch.jhu.edu) cluster using the login credentials we provided on monday. It will not complete on your local machine.__
# 
# We will use slurm to run the read qc program `fastqc` on a collection of fastq files provided for you.
#
# ## Learning Objectives
# - Familiarize yourself with the basic slurm commands
# - Learn how to submit jobs to the cluster using the slurm command `sbatch`
# - Learn how to monitor jobs using the slurm command `sacct`
# - Understand where and how to find the output of your jobs.
# - Gain experience in reviewing the output of fastqc to assess the quality of the reads in a fastq file.
#

#%% [markdown]
# ## The Task
# We have provided you with a set of raw read fastq files in the directory `/data/me440_lgoff2/datasets/` on the rockfish cluster.  Your task is to copy the fastq files to your working directory and run `fastqc` on each of these files and review the output to assess the quality of the reads in each file.  
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
# The samples for this study were each sequenced once (1 run per sample) on the Illumina sequencing platform to generate single-end reads of 100bp in length.
# 
# The directory that contains the gzip-compressed .fastq.gz files is `/data/me440_lgoff2/datasets/RNA-Seq/data/raw/GSE74985`. 
#
# Let's use the `cp` command to copy all of the .fastq.gz files from the shared directory to our new local directory.

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
# For this exercise we want to run the read quality control program `fastqc` on each of these files.  
# 
# `fastqc` is a program that parses the reads in a fastq file and outputs some summary statistics, metrics, and heuristics to tell us more about what is in the file, and the quality of the reads.  It generates a report that can be used to identify systematic issues/errors that can be common in the library preparation or sequencing of the reads.  
#
# Let's run fastqc to see if we can get an idea of the quality of the reads in one of the files. 
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
# The `-o` flag tells `fastqc` where to store the output files, in this case, the directory `fastqc_output`.
#
# _Remember, for most tools we can either use the `man` command or the `--help` flag to get more information about the tool and how to use it. `fastqc` is no exception here and if we want to find other arguments that are available we can run:_
#
# `fastqc --help`
# Let's try the first file `SRR2916027.fastq.gz`

#%%
%%bash
fastqc -o fastqc_output data/SRR2916027.fastq.gz

#%% [markdown]
# This generates a report for the file `SRR2916027.fastq.gz` and stores it in the directory `fastqc_output` and takes ~3 minutes to run.
# 
# Let's take a look at the .html report that was generated. To do so, we can either download the file and open in a browser, or preview the file in vscode.  Let's try the latter.
#
# (Cmd+Shift+P) -> HTML: Open Preview (_You may have to install the `HTML Preview` extension for VSCode_)

#%% [markdown]
# We still have a number of files for which to run this report.  We could run each job one at a time, but that would take a lot of babysitting.
#
# We could _also_ write a for loop to run each job in succession but it would still run each job sequentially (~3min * 24 jobs = ~72 minutes).
#
# Instead, we can use slurm to submit each job to the cluster and let slurm manage the resources, scheduling, and execution of each job in parallel.

#%% [markdown]
### SLURM: An Introduction
#
# [Slurm](https://slurm.schedmd.com/quickstart.html) (Simple Linux Utility for Resource Management) is an open-source job scheduler that allocates compute resources on clusters for queued user jobs. 
# Slurm has become a standard for supercomputing environments, providing both resource management and job scheduling. Slurm is used on the [rockfish cluster](http://edulogin.arch.jhu.edu).
# 
# Slurm is a command-line tool that can be used to submit, monitor, and cancel jobs. 
#
# Slurm is _primarily_ useful when we need to run a large number of the same type of jobs that can be run independently of each other. This type of problem is called 'embarrasingly parallel'. Instead of running each of the jobs in succession, we can submit them all at once and let Slurm manage the resources, scheduling, and execution of each job in parallel.
# 
# This is the case for many bioinformatics pipelines, including for example, the alignment of many samples to a reference genome. 

#%% [markdown]
#### Submitting Jobs with `sbatch`
#
# `sbatch` is used to submit a job to the scheduler for execution. 
# A 'job' is a set of commands that you would like to execute. 
#
# Typically, you provide `sbatch` with a script that provides directives and commands, but for simpler use cases like this, the `--wrap` argument allows for direct submission of command-line calls.
#
# #### Submission using `--wrap`:
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
# Notice the only information that we get from this is an acknowledgement that the job was submitted and a `job id` number.
#
# By default, the output of the job is written to a file called `slurm-<job_id>.out` in the current working directory. This file captures the STDOUT produced by the job (what is normally printed to the terminal.) 
#
# Because the job was submitted and executed on the cluster, the output file is written to this file so we have a record of what _would_ have been printed to the screen, in this case, the string 'Hello World!'.
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
# The job has now been 'submitted' to the slurm scheduler and we have a new `job_id` number. We can open the corresponding `.out` file to watch as the output that was previously printed to the screen is written to the file. 
#
# It's also important to note that this no longer ties up our available resources in the terminal to run this job.  It's actually being executed somewhere else on the cluster.  So we can keep working while the job is being executed.
#

#%% [markdown]
# We can check the status of our jobs using the slurm command `sacct` ([Slurm 'Accounting'](https://slurm.schedmd.com/sacct.html)).  Let's take a look.

#%%
%%bash
sacct

#%% [markdown]
# We can see that our job is currently running.  
# 
# We can also see the job id number, the user who submitted the job, the start time, the partition, the state, and any exit codes (errors) that might have been produced by the job.
#

#%% [markdown]
# Ok, let's try and save ourselves some time and parallelize the submition and execution of each of the fastqc jobs for the `.fastq.gz` files in our `data` directory.
#
# We could write out each of the jobs that we want to execute and wrap each in an `sbatch` call, but that could be tedious and error prone.
# 
# Instead, we can use a bash `for` loop to iterate over each of the files in our `data` directory and submit a job for each file.
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
# Now we should see a lot more jobs in the queue. We can see that some of the jobs are running, some may have completed, and some may still be 'PENDING'.
# 
# This is because slurm is 'balancing' the available resources on the cluster to ensure that all jobs are able to complete in a timely manner.  This means that some jobs may be queued until there are enough resources available to execute them.
#
# In general, slurm is a 'first-come, first-served' model in which jobs that are submitted first will have more priority for execution than those submitted later. 
#

#%% [markdown]
# ### Monitoring All Jobs with `squeue`
#
# While `sacct` can be useful for monitoring your own jobs, `squeue` displays the status of all 'active' jobs on the cluster, including information about those submitted from other users.
# 

#%% [bash]
%%bash
squeue

#%% [markdown]
# We can see that there are a lot of jobs in the queue.  We can also see that the jobs are being distributed across the different compute nodes (NODELIST) on the cluster.
#
# This can give you a feel for how busy the cluster is and how many resources might be available for your jobs.
#
# To focus again on _your_ specific jobs, we can use the `-u` flag to filter by user.

#%%
%%bash
squeue -u lgoff2

#%% [markdown]
#### Cancelling Jobs with `scancel`
#
#`scancel` is used to terminate 'RUNNING' or 'PENDING' jobs.
#
##### Cancel a Specific Job:
# You can cancel a specific job by providing the job id number to the `scancel` command:
#
# `scancel <job_id>`

#%% [markdown]
##### Cancel All Jobs for a User:
# Additionally, if you found you made a mistake while submitting a large number of jobs, you can cancel all jobs for a specific user:
#
# `scancel -u lgoff2`
# 
# _A word of caution however, this will cancel **all** of your jobs, including the job that was submitted when we created our ssh tunnel to allow VSCode to connect.  Doing this will boot you off of the cluster and you will have to reconnect._

#%% [markdown]
#### Inspecting Cluster with `sinfo`
# As we learn more about how to use the cluster and slurm, it may be useful to learn more about how the cluster is configured and what resources are available for slurm.
#
# `sinfo` provides an overview of Slurm nodes (compute nodes) and partitions ('queues' for job submission).
# 

#%% 
%%bash
sinfo

#%% [markdown]
# This tells us that there is one (unique) partition or 'queue' for job submission called `defq` (for 'default queue'). It is currently available (up) with a 2 HR timelimit per job (`2:00:00`).
#
# It also tells us the number of compute nodes available for this partition (6) and the 'NODELIST' gives us the 'names' of the compute nodes that are included as part of this partition. 


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
#### Running Jobs Interactively with `srun`
# We've already seen `srun` in action when we used it to create an interactive session and the ssh tunnel on the cluster.
#
# Let's revisit this command again to deconstruct what it is doing.
#
# `ssh lgoff2@edulogin.arch.jhu.edu "srun --time=2:00:00 --mem-per-cpu=4GB --cpus-per-task=2 VSCode-linux-x64/bin/code tunnel --accept-server-license-terms"`
#
# The above command was used to `ssh` into the cluster with my credentials (lgoff2@edulogin.arch.jhu.edu), and then immediately excute a call to `srun` on the login node.
#
# `srun` was used to initiate an interactive session on the cluster with the following resources:
# - `--time=2:00:00` - 2 hours of walltime
# - `--mem-per-cpu=4GB` - 4 Gb of memory per cpu
# - `--cpus-per-task=2` - 2 cpus per task
#
# By default `srun` will start an interactive session and run whatever commands are provided.  He're asking `srun` to create an interactive session and _immediately_ run the following command in the interactive session:
#
# `VSCode-linux-x64/bin/code tunnel --accept-server-license-terms`
#
# This will start the VSCode server on the cluster and create the ssh tunnel which allows us to connect to the server.
#


#%% [markdown]
# Often times when we're debugging and or fleshing out job submissions, having an interactive session on the cluster is helpful.
#

#%% [markdown]
# ## MultiQC
# Ok, back to the matter at hand, we have successfully executed `fastq` on each of our fastq files and generated a report for each file.
# 
# We'd like to actually start to look through each of these reports to see if anything looks odd, or otherwise suspect in terms of read quality for these samples.
#
# But parsing through 24 .html files to find patterns/trends might be a bit tedious.
#
# [MultiQC](https://multiqc.info/) is an exceptionally useful tool for aggregating the output of multiple steps or samples in a bioinformatics workflow (including fastqc runs) into a single report.
#
# MultiQC traverses a directory and searches for output from common bioinformatics tools, extracts the data, and generates a single report summarizing the data.
#
# Let's use MultiQC to aggregate the output of our fastqc runs into a single report that will be much easier to make comparisons across our samples.

#%%
%%bash
multiqc .

#%% [markdown]
# Let's take a look at the summary report that was generated.
#
# See anything useful or interesting about the samples?

# %% [markdown]
#  # Module 5 Problem Set
#  Name: **<Your Name>**

# %% [markdown]
#  For this week's problemset, we will be looking at the following study:
# 
#  * Dumitriu, A., Golji, J., Labadorf, A.T. et al. Integrative analyses of proteomics and RNA transcriptomics implicate mitochondrial processes, protein folding pathways and GWAS loci in Parkinson disease. BMC Med Genomics 9, 5 (2015). [https://doi.org/10.1186/s12920-016-0164-y](https://doi.org/10.1186/s12920-016-0164-y)
# 
#  Summary:
# 
#  Parkinson disease (PD) is marked by the accumulation of specific proteins in neuron aggregates, with previous studies focusing mainly on DNA sequence variants and RNA levels.
#  In this research, both state-of-the-art proteomics and RNA-sequencing technologies were employed to analyze prefrontal cortex samples from PD patients and healthy controls.
#  The study identified significant differences in proteins and genes between PD and controls, revealing insights into mitochondrial-related pathways and protein folding pathways, and highlighted the benefit of combining multiple genome-wide platforms to gain a deeper understanding of PD pathology.
# 
#  You will be locating the data and metadata associated with this study, and performing some basic quality control analyses on a subset of the RNA-Seq data.
# 
# 

# %% [markdown]
#  **What is the GEO accession number for this study?**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  **Briefly describe the different types of data that the authors have made available under this study accession number?**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
# ### Sample Metadata
#  Visit the GEO (Gene Expression Omnibus) database and locate the GEO series with the above accession number.
# 
#  **What is the url of the GEO series for this dataset?**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  **How many samples (in total, all samples) are present in this GEO series?**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  Open the 'SRA Run Selector' link for this series and download the metadata (sample information) for this study.
# 
#  **Load the metadata (the text file is actually a .csv) into a pandas dataframe**, show your work below and **print the first 5 rows of the dataframe**.

# %%
import pandas as pd

# Your code here:


# %% [markdown]
#  We are going to focus on one specific sample. The sample we are after has a `BioSample` id of `SAMN03651528`.
# 
#  `BioSample` is one of the columns in the metadata dataframe you just created.
# 
#  Fiter your dataframe to find the `Run` column value (SRR accession number) for the sample with the `BioSample` id of `SAMN03651528`.
# 

# %%
# Your code here:



# %% [markdown]
# ### Fetching a Fastq File
#  Next you will locate and download the .fastq.gz file associated with the `BioSample` named `SAMN03651528`.
# 
#  We have provided a function below that will download the file for you, given the SRA `RUN` accession number (SRR...) for a given sample.
# 
#  You will need to provide the SRA `RUN` accession number for the sample you selected above to the function below.
# 
#  **Note:** This file is ~2.2Gb so it may take some time to download and will take up space on your computer.

# %%
from urllib.request import urlretrieve

def download_fastq(SRRAccession):
    baseURL = "https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc="
    fname = SRRAccession + ".fastq.gz"
    path, headers = urlretrieve(baseURL + SRRAccession, fname)
    return path

# Use the function above to download the fastq.gz file for the 'SAMN03651528' sample using it's SRR accession number.

SRR_ID = ''
download_fastq(SRR_ID)


# %% [markdown]
# ### FastQC Analysis
#  After downloading the FASTQ file, perform a quality control check using FastQC.
# 
#  If you need to install fastqc on your computer, you can do so using mamba:
# 
#  `mamba install -c bioconda fastqc`

# %%
# Run fastqc on the fastq file you downloaded above.
%%bash


# %% [markdown]
#  Use the information in the fastqc .html report to answer the following questions
# 

# %% [markdown]
#  **How many Total Sequences are in the fastq file for this sample?**

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  **What is the sequence length of the reads?**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  **Which two metrics are indicating warnings (yellow)? For each, describe what you think may be causing the deviation from the expected result.**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  **Which two metrics are indicating significant issues (red)?**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  Take a look at the "Overrepresented sequences" section.
# 
#  **List any overrepresented sequences.**
# 

# %% [markdown]
# 
# 

# %% [markdown]
#  **Describe how these overrepresented sequences may have arisen.**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  Both the 'Adapter Content' plot, and the 'Per base sequence content' plot indicate that there may be deviation from expected sequence content.
# 
#  In fact they are both indicating the same issue.
# 
#  **What is the issue and what do you think may be causing it?**
# 

# %% [markdown]
#  Answer:
# 

# %% [markdown]
#  **How might you go about fixing this issue for the affected reads?**
# 

# %% [markdown]
#  Answer:
# 



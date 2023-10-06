#%% [markdown]
# # Module 5 Problem Set
# Name: **_Your Name_**

#%% [markdown]
# For this week's problemset, we will be looking at the following study:
#
# * Dumitriu, A., Golji, J., Labadorf, A.T. et al. Integrative analyses of proteomics and RNA transcriptomics implicate mitochondrial processes, protein folding pathways and GWAS loci in Parkinson disease. BMC Med Genomics 9, 5 (2015). [https://doi.org/10.1186/s12920-016-0164-y](https://doi.org/10.1186/s12920-016-0164-y)
#
# Summary:
#
# Parkinson disease (PD) is marked by the accumulation of specific proteins in neuron aggregates, with previous studies focusing mainly on DNA sequence variants and RNA levels.
# In this research, both state-of-the-art proteomics and RNA-sequencing technologies were employed to analyze prefrontal cortex samples from PD patients and healthy controls. 
# The study identified significant differences in proteins and genes between PD and controls, revealing insights into mitochondrial-related pathways and protein folding pathways, and highlighted the benefit of combining multiple genome-wide platforms to gain a deeper understanding of PD pathology.
#
# You will be locating the data and metadata associated with this study, and performing some basic quality control analyses on a subset of the RNA-Seq data.
#
#
# ## 1. Data Access:
#
# 

#%% [markdown]
# 1.1 What is the GEO accession number for this study?  Describe the different types of data that the authors have made available under this study accession number?

#%% [markdown]
#
#

#%% [markdown]
# Visit the GEO(Gene Expression Omnibus) database and locate the GEO series with the above accession number.  What is the url of the GEO series for this dataset?

#%% [markdown]
# How many samples (in total, all samples) are present in this GEO series?

#%% [markdown]
#

#%% [markdown]
# Open the 'SRA Run Selector' for this series and download the metadata (sample information) for this study. 
#
# Load the metadata into a pandas dataframe, show your work below and print the first 5 rows of the dataframe.

#%%
import pandas as pd

#%% [markdown]
#
#


#%% [markdown]
# Locate and download the .fastq file associated with the sample named `C_0003`. We have provided a function below that will download the file for you, given the SRA `RUN` accession number (SRR...) for a given sample. 
#
# You will need to provide the SRA `RUN` accession number for the sample you selected above to the function below.
# 
# **Note:** This file is ~2.2Gb so it may take some time to download and will take up space on your computer.

#%%
from urllib.request import urlretrieve

def download_fastq(SRRAccession):
    baseURL = "https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc="
    fname = SRRAccession + ".fastq.gz"
    path, headers = urlretrieve(baseURL + SRRAccession, fname)
    return path

# Use the function above to download the fastq.gz file for the sample you selected

    
#%% [markdown]
  
    

#%% [markdown]
# 2. FastQC Analysis
# After downloading the FASTQ file, perform a quality control check using FastQC.
# 
# a. Provide a summary of the basic statistics tab.
# b. Which regions of the reads, if any, are flagged as having poor quality?
# c. Interpret the "Per base sequence content" graph. Are there any irregularities? What could they indicate?
# d. Comment on the "Overrepresented sequences" section. Are there any overrepresented sequences? If so, what might be the cause?

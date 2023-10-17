#%% [markdown]
# ## Install pydeseq2

#%% [bash]
%%bash
conda install -c bioconda pydeseq2

# %%
# import pydeseq2 
import pandas as pd
import numpy as np

#%% [markdown]
# PyDESeq2 requires two types of inputs:
#
# - A count matrix of shape ‘number of samples’ x ‘number of genes’, containing read counts (non-negative integers),
#
# - Metadata (or “column” data) of shape ‘number of samples’ x ‘number of variables’, containing sample annotations that will be used to split the data in cohorts.

#%% 

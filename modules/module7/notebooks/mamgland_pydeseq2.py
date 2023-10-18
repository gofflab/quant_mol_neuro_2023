#%% [markdown]
# # Module 7 - Differential Expression Analysis of bulk RNA-Seq
# In this notebook, we will be showing and example differential expression analysis using the pydeseq2 package.
#
# `pydeseq2` is a python implementation for the popular R package `DESeq2` for differential analysis of **raw** count data from bulk RNA-Seq experiments.
#
# We will be using the following tools/packages:
#
# - [`pydeseq2`](https://pydeseq2.readthedocs.io/en/latest/) for differential expression analysis
# - [`anndata`](https://anndata.readthedocs.io/en/latest/) for data wrangling and to generate the base object required by `pydeseq2`
# - [`scanpy`](https://scanpy.readthedocs.io/en/latest/) for select visualizations of the pydesq2 object
# - [`plotnine`](https://plotnine.readthedocs.io/en/latest/) for other visualizations
# - [`seaborn`](https://seaborn.readthedocs.io/en/latest/) for heatmap (clustermap) visualizations
# - `pandas` and `numpy` for data wrangling
#
#  Much of the code in this notebook is adapted from the official [pydeseq2 'Step-by-step' workflow](https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_step_by_step.html#sphx-glr-auto-examples-plot-step-by-step-py).
# 
# First, let's install the required packages
#
# ## Install required packages
# You will only need to run this cell once!
#
# Here we are also using pip to install the most recent version of scanpy from github.

#%% [bash]
%%bash
#mamba install -c bioconda pydeseq2 leidenalg seaborn anndata
#pip install git+https://github.com/scverse/scanpy

# %% [markdown]
# Now we can import all of the required packages

# %%
import pandas as pd
import numpy as np
import plotnine as pn
import seaborn as sns

import anndata as ad

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scanpy as sc

#%% [markdown]
# ## Load Data
# The three required input files for this analysis are provided in the `data` directory.
# 
# PyDESeq2 can use either the raw count matrix and the metadata table, or (more conveniently) an `anndata` object containing all of the feature (`adata.var`) and sample information (`adata.obs`) as well as the count matrix (`adata.X`).
#
# We are going to load each of the three required matrices using `pandas.read_csv`

#%% 
counts = pd.read_csv("data/mamgland_counts.txt",index_col=0)
sampleData = pd.read_csv("data/mamgland_colData.txt",index_col=0)
geneData = pd.read_csv("data/mamgland_rowData.txt",index_col=0)

#%% [markdown]
# Let's take a quick look at each of the matrices to make sure they were imported and formatted correctly.

# %%
counts.head()

# %%
sampleData.head()

# %%
geneData.head()

# %% [markdown]
# Because the counts matrix is being provided as a `genesXsamples` matrix, we need to transpose it to a `samplesXgenes` matrix to create an `anndata` object.
#

# %%
# Transpose the count matrix to samples x genes
counts = counts.T

# %% [markdown]
# Now that we have the required matrices, we can create the `anndata` object `adata`

#%% Create AnnData Object
adata = ad.AnnData(X=counts,obs=sampleData,var=geneData)

#%%
adata

#%% [markdown]
# ## Read counts modeling
# Read counts modeling with the DeseqDataSet class
#
#The DeseqDataSet class takes one required positional argument (the `anndata` object), as well as a set of optional keyword arguments, including:
#
# - design_factor: the name of the column of metadata to be used as a design variable
#
# - refit_cooks: whether to refit cooks outliers – this is advised, in general.
#
# We construct the DeseqDataSet object from `adata` and specify the design factor to be the `group` column of the metadata table, which contains the information about the sample groups (Basal, LP, and ML).

# %%
dds = DeseqDataSet(
    adata = adata,
    design_factors="group",  # compare samples based on the "condition"
    # column ("B" vs "A")
    refit_cooks=True,
    n_cpus=8,
)

#%%
dds

#%% [markdown]
# Notice that the `dds` object is actually an 'extended' `anndata` object, with all of the original data slots, accessors, methods, etc. The only difference is the addition of a few additional attributes and methods which will be used for the specific purpose of differential expression analysis.

#%% [markdown]
# ## Some Basic QC
# Let's take a quick look at the data to make sure everything looks good.
# One of the ways we can do this is by performing a principal component analysis (PCA) and visualizing the results.
#
# PCA is a dimensionality reduction technique that is often used to visualize high-dimensional data. And can be useful to identify relationships between samples and/or groups of samples as well as identify potential outliers.
#
# We can use the `scanpy` package to perform the PCA and visualize the results. `scanpy` is a python package that is originally designed for single-cell RNA-Seq analysis, but can adopt some of the qc and visualization tools for bulk RNA-Seq data as well.
#
#%% PCA
sc.tl.pca(dds)
sc.pl.pca(dds, color="group", title="PCA", s=100)
# %%

#%% [markdown]
# ### Compute normalization factors

#%%
dds.fit_size_factors()
dds.obsm["size_factors"]

#%% [markdown]
# ### Fit gene-wise dispersions

#%%
dds.fit_genewise_dispersions()

dds.varm["genewise_dispersions"]

# %% [markdown]
# ### Fit dispersion trend coefficients

#%%
dds.fit_dispersion_trend()
dds.uns["trend_coeffs"]
dds.varm["fitted_dispersions"]

# %% [markdown]
# ### Dispersion priors

#%% 
dds.fit_dispersion_prior()
print(
    f"logres_prior={dds.uns['_squared_logres']}, sigma_prior={dds.uns['prior_disp_var']}"
)

# %% [markdown]
# ### MAP Dispersions
# The fit_MAP_dispersions method filters the genes for which dispersion shrinkage is applied. 
# 
# Indeed, for genes whose MLE dispersions are too high above the trend curve, the original MLE value is kept. 
# 
# The final values of the dispersions that are used for downstream analysis is stored in dds.dispersions.

#%%
dds.fit_MAP_dispersions()
dds.varm["MAP_dispersions"]
dds.varm["dispersions"]

# %% [markdown]
# ### Fit log fold changes
# Note that in the DeseqDataSet object, the log-fold changes are stored in natural log scale, but that the results dataframe output by the summary method of DeseqStats displays LFCs in log2 scale (see later on).

#%%
dds.fit_LFC()
dds.varm["LFC"]

# %% [markdown]
# ### Calculate Cooks distances and refit
# Note that this step is optional.

#%%
dds.calculate_cooks()
if dds.refit_cooks:
    # Replace outlier counts
    dds.refit()

# %% [markdown]
# ## Statistical analysis
# Statistical analysis with the DeseqStats class. The DeseqDataSet class has a unique mandatory arguments, dds, which should be a fitted DeseqDataSet object, as well as a set of optional keyword arguments, among which:
#
# - alpha: the p-value and adjusted p-value significance threshold
#
# - cooks_filter: whether to filter p-values based on cooks outliers
#
# - independent_filter: whether to perform independent filtering to correct p-value trends.
# %%
stat_res = DeseqStats(dds, alpha=0.05, cooks_filter=True, independent_filter=True)

# %% [markdown]
# ### Wald tests
# 

#%%
stat_res.run_wald_test()
stat_res.p_values
# %% [markdown]
# ### Cook's filtering
# Note: this step is optional

#%% 
if stat_res.cooks_filter:
    stat_res._cooks_filtering()
stat_res.p_values

# %% [markdown]
# ### P-value adjustment

#%%
if stat_res.independent_filter:
    stat_res._independent_filtering()
else:
    stat_res._p_value_adjustment()

stat_res.padj

# %% [markdown]
# ### Building a results dataframe
# This dataframe is stored in the results_df attribute of the DeseqStats class.

#%%
stat_res.summary()

# %% [markdown]
# Let's store the summary output into it's own dataframe for easier access.

#%%
res = stat_res.results_df.copy()

# %% [markdown]
# ### LFC Shrinkage
# For visualization or post-processing purposes, it might be suitable to perform LFC shrinkage. This is implemented by the lfc_shrink method.

#%%
stat_res.lfc_shrink(coeff="group_LP_vs_Basal")

#%% [markdown]
# ## Visualizations

#%%
plot_data = stat_res.results_df.copy()
plot_data["neg_log_10_padj"] = -np.log10(plot_data["padj"])
pval_thres = 1e-20
plot_data["significant"] = "Not significant"
plot_data["significant"][plot_data["padj"] <= pval_thres] = "Significant"

#%%
volcano_plot = (
    pn.ggplot(plot_data,pn.aes(x="log2FoldChange",y="neg_log_10_padj"))
    + pn.geom_point(pn.aes(color="significant"),alpha=0.5)
    + pn.scale_color_manual(values=["black","red"])
)

volcano_plot

# %% [markdown]
# ### Heatmap of significant genes

#%%
dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])

# %%
sig_gene_df = res[(res['padj'] <= pval_thres) & (abs(res["log2FoldChange"] > 0.5)) & (res["baseMean"] > 20 )]

#%%
dds_sig_genes = dds[:,sig_gene_df.index]

# %%
plot_data = pd.DataFrame(dds_sig_genes.layers['log1p'].T, index = dds_sig_genes.var["SYMBOL"], columns = dds_sig_genes.obs_names)

# %%
sns.clustermap(plot_data, z_score=0,cmap="RdYlBu_r")

# %%

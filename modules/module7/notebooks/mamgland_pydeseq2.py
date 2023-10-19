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
# Subset `adata` to only include the samples we want to compare (['Basal','LP'])

#%%
adata = adata[adata.obs['group'].isin(['Basal','LP'])].copy()

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

#%% [markdown]
# Next we can check to make sure that the variance is stabilized appropriately.
#
# The `plot_dispersions` method plots the estimated, shrunken, and fitted dispersion values against the mean of normalized counts.

#%%
dds.plot_dispersions()

# %% [markdown]
# ### Fit log fold changes
# Note that in the DeseqDataSet object, the log-fold changes are stored in natural log scale, but that the results dataframe output by the summary method of DeseqStats displays LFCs in log2 scale (see later on).
#

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
# Statistical analysis with the DeseqStats class. We will create an instance of `DeseqStats` called `stat_res` and pass the `dds` object to it.
#
# This object is responsible for performing the statistical analysis and storing the results.
#
# In addition to the required `dds` input, there are a few additional optional keyword arguments, including:
#
# - alpha: the p-value and adjusted p-value significance threshold
#
# - cooks_filter: whether to filter p-values based on cooks outliers
#
# - independent_filter: whether to perform independent filtering to correct p-value trends.
# %%
stat_res = DeseqStats(dds, alpha=0.05, cooks_filter=True, independent_filter=True)

# %% [markdown]
# Let's inspect the `stat_res` object to see what it contains.

#%%
stat_res.design_matrix

#%%
stat_res.contrast

#%%
stat_res.contrast_vector

#%%
stat_res.LFC

#%% [markdown]
# A few additional slots will be filled after the statistical analysis is performed, including:
# - p_values: the p-values of the Wald tests
# - padj: the adjusted p-values
# - statistics: the Wald test statistics
#
# And a few others.  For now, let's use `stat_res` to actually perform the statistical analysis.

# %% [markdown]
# ### Wald tests
# In bulk RNA-Seq with the DESeq2 workflow, the Wald test evaluates the probability of the null hypothesis that the log fold change between two conditions is equal to 0.
#
# In `pydeseq2`, the Wald test is performed by the `run_wald_test` method of the `DeseqStats` class (in our case, the instance `stat_res`).

#%%
stat_res.run_wald_test()
stat_res.p_values

# %% [markdown]
# ### Cook's filtering
# Note: this step is optional
# 
# Cook's filtering is used in this context to remove genes with extreme counts that may be driving the results.
#
# The `cooks_filter` method of the `DeseqStats` class performs Cook's filtering by removing genes with Cook's distance above a certain threshold.

#%% 
if stat_res.cooks_filter:
    stat_res._cooks_filtering()

# %% [markdown]
# Now that we've performed the Wald test and Cook's filtering, we actually have raw p-values associated with the Wald test statistics for each gene.

#%%
stat_res.p_values

# %% [markdown]
# ### P-value adjustment
# Because we have performed a large number of statistical tests (one for each gene), we need to adjust the p-values to account for multiple testing.
# 
# The `_p_value_adjustment` method of the `DeseqStats` class performs the p-value adjustment.

#%%
stat_res._p_value_adjustment()

#%% [markdown]
# The results of the p-value adjustment are stored in the `padj` attribute of the `stat_res` object.
# 
#%%
stat_res.padj

# %% [markdown]
# ### Building a results dataframe
# To put this all together into something more easily interpretable, we can build a results dataframe that contains the results of our statistical analysis.
#
# The results_df attribute of the DeseqStats class can store this dataframe once it's been generated. 
# 
# To generate this results dataframe, we have to first make a call to the `summary()` method of the DeseqStats class.

#%%
stat_res.summary()

# %% [markdown]
# Let's store the summary output into it's own dataframe for easier access.

#%%
res = stat_res.results_df.copy()

#%% [markdown]
# Now let's take a look a the p-value distribution to see if this looks appropriate

#%%
pval_hist = (
    pn.ggplot(stat_res.results_df, pn.aes(x="pvalue"))
    + pn.geom_histogram()
    + pn.ggtitle("P-value distribution")
)

pval_hist

# %% [markdown]
# ### LFC Shrinkage
# For visualization or post-processing purposes, it might be suitable to perform LFC shrinkage. This is implemented by the lfc_shrink method.

#%%
stat_res.lfc_shrink(coeff="group_LP_vs_Basal")

#%% [markdown]
# ## Visualizations
# Let's explore some of our differential analysis results using some visualizations.
#
# First, we will make a volcano plot of the log2 fold change vs the -log10 adjusted p-value.
#
# We will use the shrunken LFCs for the x axis and a calculated -log10 `padj` for the y axis.
#
# Below, we are also choosing a p-value threshold of 1e-20 to highlight the most significant genes. 

#%%
# Create a copy of the results dataframe to use for plotting
plot_data = stat_res.results_df.copy()

# Calculate the -log10 padj
plot_data["neg_log_10_padj"] = -np.log10(plot_data["padj"])

# Set a p-value threshold
pval_thres = 1e-20

# Add a column to indicate whether a particular gene test is significant or not
plot_data["significant"] = "Not significant"
plot_data["significant"][plot_data["padj"] <= pval_thres] = "Significant"

#%% [markdown]
# Now that we've prepared the plotting data, we can make the volcano plot using `plotnine`

#%%
volcano_plot = (
    pn.ggplot(plot_data,pn.aes(x="log2FoldChange",y="neg_log_10_padj"))
    + pn.geom_point(pn.aes(color="significant"),alpha=0.5)
    + pn.geom_vline(xintercept=0,linetype="dashed")
    + pn.geom_hline(yintercept=-np.log10(pval_thres))
    + pn.scale_color_manual(values=["black","red"])
    + pn.ggtitle("Volcano Plot of Differential Expression - (LP vs Basal)")
)

volcano_plot

# %% [markdown]
# Is there a relationship between mean expression and p-value? Let's investigate using a scatter plot.

#%%
mean_v_pval_plot = (
    pn.ggplot(plot_data,pn.aes(x="baseMean",y="neg_log_10_padj"))
    + pn.geom_point(pn.aes(color="significant"),alpha=0.5)
    + pn.scale_color_manual(values=["black","red"])
    + pn.scale_x_log10()
    + pn.ggtitle("Mean Expression vs P-value")
)

mean_v_pval_plot

# %% [markdown]
# ### Heatmap of significant genes
# Finally, let's take a subset of our significant genes and make a heatmap of the log1p normalized counts.
#
# Since we have a large number, we're going to subset using some fairly stringent criteria:
# - A very low adjusted pvalue (`padj` <= 1e-20)
# - A log fold change  > 0.5 (`log2FoldChange`)
# - Mean expression > 20 (`baseMean`)
#
# After applying these filters, we are going to identify the genes that meet these criteria, and create an index so that we can subset the original `dds` object to only include these genes.
#
# We will then use the `sns.clustermap` function to create a heatmap of the log1p normalized counts for these genes.
#
# First, lets create a new layer in `dds` that contains the log1p normalized counts.

#%%
dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])

#%% [markdown]
# Now, let's apply our filtering critera to the stat_res.results_df dataframe to identify the genes that meet our criteria.
# %%
sig_gene_df = stat_res.results_df[(stat_res.results_df['padj'] <= pval_thres) & (abs(stat_res.results_df["log2FoldChange"] > 0.5)) & (stat_res.results_df["baseMean"] > 20 )]

# %% [markdown]
# `sig_gene_df` is a new results dataframe that only contains the genes that meet our criteria.

#%%
sig_gene_df.shape

# %% [markdown]
# Now we'll use the index of `sig_gene_df` (which contains the gene identifiers) to subset the `dds` object to only include the genes that meet our significance criteria.
#%%
dds_sig_genes = dds[:,sig_gene_df.index]

#%%
dds_sig_genes

# %% [markdown]
# Now, we will take the `dds_sig_genes` object and create a new `wide` dataframe that contains the log1p normalized counts for each sample.
#
# This dataframe will have the shape (sigGenes X samples).

# %%
plot_data = pd.DataFrame(dds_sig_genes.layers['log1p'].T, index = dds_sig_genes.var["SYMBOL"], columns = dds_sig_genes.obs_names)

#%% [markdown]
# And finally, we can use the seaborn `clustermap` function to create a heatmap of the log1p normalized counts for these genes.
# 
# The argument `z_score=0` will normalize the rows (genes) to have a mean of 0 and a standard deviation of 1 to allow for easier visualization of the differences between samples.

# %%
sns.clustermap(plot_data, z_score=0,cmap="RdYlBu_r")

# %%

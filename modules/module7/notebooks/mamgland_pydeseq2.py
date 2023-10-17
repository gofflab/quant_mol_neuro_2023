#%% [markdown]
# ## Install pydeseq2

#%% [bash]
%%bash
#mamba install -c bioconda pydeseq2 leidenalg
#pip install git+https://github.com/scverse/scanpy

# %%
# import pydeseq2 
import pandas as pd
import numpy as np
import plotnine as pn

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scanpy as sc

#%% [markdown]
# PyDESeq2 requires two types of inputs:
#
# - A count matrix of shape ‘number of samples’ x ‘number of genes’, containing read counts (non-negative integers),
#
# - Metadata (or “column” data) of shape ‘number of samples’ x ‘number of variables’, containing sample annotations that will be used to split the data in cohorts.

#%% 
counts = pd.read_csv("data/mamgland_counts.txt",index_col=0)
sampleData = pd.read_csv("data/mamgland_colData.txt",index_col=0)
geneData = pd.read_csv("data/mamgland_rowData.txt",index_col=0)

# %%
counts.head()

# %%
sampleData.head()

# %%
geneData.head()

# %%
# Transpose the count matrix to samples x genes
counts = counts.T

#%% [markdown]
# ## Read counts modeling
# Read counts modeling with the DeseqDataSet class
#
#The DeseqDataSet class has two mandatory arguments, counts and metadata, as well as a set of optional keyword arguments, among which:
#
# - design_factor: the name of the column of metadata to be used as a design variable
#
# - refit_cooks: whether to refit cooks outliers – this is advised, in general.

# %%
dds = DeseqDataSet(
    counts=counts,
    metadata=sampleData,
    design_factors="group",  # compare samples based on the "condition"
    # column ("B" vs "A")
    refit_cooks=True,
    n_cpus=8,
)

#%% [markdown]
# ## Some Basic QC
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

# %%

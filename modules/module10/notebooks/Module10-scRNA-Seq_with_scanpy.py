#%% [markdown]
# # Module 10 - Single Cell RNA-Seq analysis with Scanpy
#
# ## Introduction
#Single cell analysis is a rapidly growing field in biology. The ability to analyze individual cells has allowed researchers to gain a better understanding of the heterogeneity of cells within a population. 
# 
# This module will introduce the basic concepts of single cell analysis and provide an overview of a select number of tools available for single cell analysis.

#%% [markdown]
# ## Learning Objectives
#
# - Understand the basic steps of single cell RNA-Sequencing analysis workflows
# - Develop a baseline awareness of cellular heterogeneity both between and within
#   cell 'types'.
# - Learn to identify and examine cell state transitions via pseudotime analysis
# - Understanding the application of dimensionality reduction to visualization and
#   high-dimensional sequencing data analysis.

#%% [markdown]
#
# ## Load required modules
#%%
import numpy as np
import pandas as pd
import plotnine as pn
import anndata as ad
import scanpy as sc
import scFates as scf

#%% [markdown]
# ## Data
# ### Description
# 
# The dataset we are using is
# [10x 10k neurons from an E18 mouse](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/neuron_10k_v3)
# 
# Cells for this sample are from a combined cortex, hippocampus and sub
# ventricular zone of an E18 mouse.

# - 11,869 cells detected
# - Sequenced on Illumina NovaSeq with approximately 30,000 reads per cell
# - 28bp read1 (16bp Chromium barcode and 12bp UMI), 91bp read2 (transcript), and
#   8bp I7 sample barcode

# ![10x Genomics sequencing workflow overview](https://dnacore.missouri.edu/images/10xschematic.jpg)

# Experimental questions:

# - What types of cells do we expect to find?
# - What/how many cellular 'states' do we observe at the transcriptional level?
# - How well defined are different cell types
#   - What do 'transitioning' cell types/states look like?
# - Can we identify differentially expressed and/or marker genes between cell
#   types?
# - What genes change expression over the course of neuronal development?

# ### Download dataset
#

#%% [bash]
%%bash
mkdir -p data
wget -P data/ https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_filtered_feature_bc_matrix.h5 

#%% [markdown]
# ## Overview
# The data we are using *has already been preprocessed* (aligned and demultiplexed to single cell resolution) using the Cell Ranger workflow.
#
# The `.h5` file we just downloaded is one of the output files from the 10x Genomics Cell Ranger pipeline. This file contains the following information:
#
# - A sparse matrix of gene expression values
# - A list of gene names
# - A list of cell barcodes
# - A list of cell metadata (e.g. number of genes detected, number of UMIs detected, etc.)
#   
# 
# 
# We will be using the [Scanpy](https://scanpy.readthedocs.io/en/stable/) package to analyze this data. Scanpy is a Python package for single-cell analysis. 
# 
# It has a number of useful functions for preprocessing, analyzing, and visualizing single cell data. It also has a number of useful functions for integrating single cell data from multiple samples.

# ## Set some scanpy defaults for feedback and logging
# Some default settings

#%%
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%% [markdown]
# Let's read in the datafile to create an anndata object. scanpy provides a function to directly read in 10x data from the .h5 file.
#%%
adata = sc.read_10x_h5('data/neuron_10k_v3_filtered_feature_bc_matrix.h5',gex_only=False)

#%% [markdown]
# Since we're using gene names as the index, and these are not _guaranteed_ to be unique, we'll make them unique here as a first step.

#%%
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

#%% [markdown]
# Let's take a look at the contents of the anndata object we just created.
#%%
adata

#%% [markdown]
# To access the data matrix, variable information (genes), and observation metadata (cells), we can use the `.X`, `.var`, and `.obs` attributes of the anndata object.
#
# Try them below to get a feel for what they contain.

#%%
adata.X

#%% [markdown]
# Let's start by taking a peek at some of the most highly expressed genes in these cells.
#%%
sc.pl.highest_expr_genes(adata, n_top=20)

#%% [markdown]
### Basic filtering
# Here we are going to assess some basic quality control metrics and decide which cells and genes to keep for downstream analysis.
#
# We can start by removing cells that have a low number of genes detected and genes that are detected in a low number of cells.
#
# _Note: these thresholds are relatively arbitrary and should be adjusted based on the dataset and experimental design._
#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#%% [markdown]
# Let's take a look at the data after filtering.

#%%
adata

#%% [markdown]
# Next, we will calculate some QC metrics for each cell and add them to the `.obs` attribute of the anndata object.
# 
# As part of this process, we want to assess the relative abundance of mitochondrial genes in each cell, as this can be an indicator of cell health/quality.
#
# We can do this by calculating the percentage of counts that come from mitochondrial genes for each cell.
#
# To do this, we first need to identify which genes are mitochondrial genes. We can do this by looking at the gene names and identifying which ones start with `mt-`.
#%%
# Identify mitochondrial genes by gene name
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'

adata.var

#%% [markdown]
# Now that we have identified the mitochondrial genes, we can calculate the percentage of counts that come from mitochondrial genes for each cell, as well as the total number of counts per cell and a few other metrics.
#%%
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#%% [markdown]
# This function adds a number of columns to the `.var` and `.obs` attributes of the anndata object. Let's see what's been added.
# %%
adata.var

#%% 
adata.obs

#%% [markdown]
# Let's plot a few of these metrics to see how the data look.
#
# Let's start by making a histogram of the total counts per cell and one of the number of genes detected per cell.
# 
# We can do this using the `plotnine` package.

#%%
p = (
    pn.ggplot(adata.obs,pn.aes(x='total_counts'))
    + pn.geom_histogram()
)
p.draw()

#%%
p = (
    pn.ggplot(adata.obs,pn.aes(x='n_genes_by_counts'))
    + pn.geom_histogram()
)
p.draw()

#%% [markdown]
# Here, `n_genes_by_counts` is the number of genes detected in each cell (normalized by the total counts in each cell).

#%% [markdown]
# `scanpy` also has a number of plotting functions that can be used to visualize aspects of the anndata object as well. Please be sure to look at the scanpy api for a list of available plotting functions.
#
#  Here, we're going to use the `sc.pl.violin` function to visualize the distribution of a few of the QC metrics we just calculated.
# %%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

#%% [markdown]
# And the `sc.pl.scatter` function to visualize the relationship between the number of genes detected and the total counts for each cell.

# %%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

#%% [markdown]
# Now that we have an idea of how some of the common QC metrics look, we can pick some thresholds to filter the dataset further to remove low quality cells.
#
# Let's remove cells with a low number of genes detected (<500) and a high percentage of counts that come from mitochondrial genes (>20%).
# %%
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata = adata[adata.obs.n_genes_by_counts > 500, :]

#%% [markdown]
# Redrawing these plots shows that we have removed some of the cells with low quality metrics.

# %%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# %% [markdown]
# ## Normalization
# Now that we've cleaned up the data a bit, we can move on to normalization.
#
# Normalization is the process of adjusting the data to account for differences in sequencing depth between cells and other technical effects.
#
# There are a number of different normalization methods available. Here, we will use the `normalize_total` function from `scanpy` to normalize the data by the total counts per cell.
#
# This function normalizes the counts for each cell by the total counts in that cell, multiplies the result by a scale factor (default is 10,000).

#%%
# Total-count normalize (library size correct) the data matrix
sc.pp.normalize_total(adata, target_sum=1e4)

#%% [markdown]
# After normalization, we can log transform the data to reduce the effects of high count outliers, control variance, and make the data more normally distributed.
#%%
# Log transform
sc.pp.log1p(adata)

#%% [markdown]
# ## Highly variable gene selection
# Next, we want to focus our downstream analysis on genes that are highly variable across cells, or more specifically, genes with variance in excess of what is expected based on the technical noise in the data.
#
# To do this, we fit a mean-variance relationship to the data and use this to identify genes that are outside of the expected range.

# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=4, min_disp=0.25)

#%% [markdown]
# The selection criteria above are relatively arbitrary and should be adjusted based on the dataset and experimental design.
#
# Let's take a look at the genes that were selected as highly variable.

# %%
sc.pl.highly_variable_genes(adata,)

#%% [markdown]
# Let's save the current state of the anndata object to the `raw` attribute of the anndata object. This will allow us to easily access these values later if we need to after additional filtering or processing.

# %%
adata.raw = adata

#%% [markdown]
# For the rest of the downstream analysis, we are going to only focus on the highly variable genes.
#
# So we will further subset the adata object to only include these. This will make many of the downstream steps go significantly faster (fewer genes over which to operate).

# %%
adata = adata[:, adata.var.highly_variable]

#%% [markdown]
# ## Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
# To 'correct' or 'fix' any bias that may be introduced by differences in sequencing depth between cells, we can regress out the effects of the total counts per cell and the percentage of mitochondrial genes expressed.
# 
# This should always be considered 'optional' and should be done with caution. It is not always appropriate to do this, and it is not always appropriate to do this at this stage in the analysis.
#
# For example, if you're downstream analysis requires fitting the data to a linear model, you may not want to 'fix' the data in this way, but rather, 'fit' these parameters as part of your model.  It is not appropriate to do both.
# %%
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#%% [markdown]
# ## Scale data to unit variance
# Finally, we want to scale the data to unit variance.

#%%
sc.pp.scale(adata, max_value=10)

# %% [markdown]
# ## Principal component analysis
# The next step in our workflow is to perform a dimensionality reduction on the data.
#
# Most workflows default to using principal component analysis (PCA) for this step. But other approaches may be appropriate depending on the data and the experimental design.
#
# `scanpy` provides a convenient wrapper function in it's `tools` module (`tl`) for many of these tasks including pca.

#%%
sc.tl.pca(adata, svd_solver='arpack')

#%% [markdown]
# Let's take a look at the results of the PCA. Let's start by looking at the first two principal components.

#%%
sc.pl.pca(adata)

#%% [markdown]
# We can color the plot by different attributes in the `.obs` the anndata object to identify any potential batch effects or other confounding factors.

# %%
sc.pl.pca(adata, color=["n_genes","pct_counts_mt"])

#%% [markdown]
# Here, for example, many of the cells in the top right corner of the plot have a relatively low number of genes.  This could be biological (ie erythrocytes) or it could be due to some technical issue.  
# 
# We can make a decision to go back and remove these cells from the analysis or keep this in mind as we explore further.

#%% [markdown]
# We can also look at a screeplot to see how much variance is explained by each principal component.
#
# This will help us decide how many components to include in downstream analyses and further dimensionality reductions.

# %%
sc.pl.pca_variance_ratio(adata, log=True)

#%% [markdown]
# ## Computing the neighborhood graph
# The next step in our workflow is to compute a neighborhood graph of the data. This is the first step in the process of clustering the data.
#
# The neighborhood graph is a representation of the data where each cell is connected to its nearest neighbors. This is a way of representing the data in a lower dimensional space while preserving the local structure of the data.
#
# The `sc.pp.neighbors` function computes the neighborhood graph of the data. It uses the `X_pca` attribute of the anndata object as the input data. This is why we needed to run `sc.tl.pca` before running this step.
# 
# The choice of `n_neighbors` and `n_pcs` can be modified/optimized based on the dataset and experimental design.

# %%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#%% [markdown]
# The neighbors data is stored in .uns[‘neighbors’], distances and connectivities are stored in .obsp[‘distances’] and .obsp[‘connectivities’] respectively.

#%%
adata.uns['neighbors']

#%% [markdown]
# ## Clustering the neighborhood graph
# Now that we have computed the nearest neighbor graph, we can perform a 'community detection' algorithm to identify clusters of cells along the graph.
#
# The `sc.tl.leiden` function performs Leiden community detection on the neighborhood graph. It uses the `connectivities` attribute of the anndata object as the input data.
#
# This is an iterative algorithm that can be run multiple times to identify clusters at different resolutions. The `resolution` parameter controls the granularity of the clustering. Higher values lead to more clusters.

#%%
sc.tl.leiden(adata, resolution=0.3)

#%% [markdown]
# We can revisit this clustering step later after we are able to visualize the results to decide if we need to adjust the `resolution` parameter.

#%% [markdown]
# ## Partition-based graph abstraction (PAGA)
# We can next perform partition-based graph abstraction (PAGA) to coarse-grain the graph and identify the major partitions of the dataset.
#
# This is an optional step, but it can be useful for identifying the major partitions of the dataset and visualizing the relationships between them.
# 
# [Wolf et al. (2019), PAGA: Graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells.](https://doi.org/10.1186/s13059-019-1663-x)

# %%
sc.tl.paga(adata)
sc.pl.paga(adata)

#%% [markdown]
# This gives us a course idea of how the different clusters relate to each other based on their positions in the neighbors graph.
#
# We will also use these positions to initialize the UMAP embedding in the next step.
#
# The above plot is, by default, coloring each of the clusters by their assigned `leiden` cluster number.
#%%
sc.pl.paga(adata, color='leiden')

#%% [markdown]
# We can also color the plot by other attributes in the `.obs` attribute of the anndata object, or by the average gene expression for cells in each cluster.
#
# In this way we can start to get a feel for the biological differences between the clusters
#%%
sc.pl.paga(adata, color=['Pax6','Eomes','Snap25'])

#%% [markdown]
# ## Embedding the neighborhood graph
# Now that we have computed the neighborhood graph, identified clusters of cells, set initial coarse positions for each cluster in 2D, we can embed the full graph in two dimensions to visualize the results.
#
# To do this, we will use the non-linear dimensionality reduction tool called 'UMAP' (Uniform Manifold Approximation and Projection).
# 
# Again, `scanpy` has a convenient wrapper function in it's `tools` module (`tl`) for this step.
#%%
sc.tl.umap(adata, init_pos='paga')

#%% [markdown]
# We can plot the cells in the UMAP embedding space with `sc.pl.umap`. 
# 
# Let's also color the cells by their assigned `leiden` cluster number.
# %%
sc.pl.umap(adata,color='leiden', add_outline=True)

#%% [markdown]
# We can also color the plot by other attributes in the `.obs` attribute of the anndata object, or by the expression of specific genes in each cell.

#%%
# QC metrics
sc.pl.umap(adata,color=["n_genes","pct_counts_mt"], add_outline=True)

#%%
# Progenitor and cell cycle markers
sc.pl.umap(adata,color=['Pax6','Sox2','Top2a'], add_outline=True)

# %%
# Neuronal subtype markers
sc.pl.umap(adata,color=['Eomes','Gad1','Slc17a6'], add_outline=True)

# %%
#Cortical Neuronal subtype markers
sc.pl.umap(adata,color=['Cux2','Rorb','Fezf2'], add_outline=True)

# %%
# Inhibitory neuron markers
sc.pl.umap(adata,color=['Gad1','Gad2','Slc32a1'], add_outline=True)

#%% [markdown]
# ## Finding marker genes (Differential expression w.r.t cluster)
# Now that we have identified clusters of cells, we can identify genes that are differentially expressed between clusters.
# 
# We can do this using the `sc.tl.rank_genes_groups` function. This function performs a statistical test to identify genes that are differentially expressed between groupings defined in `.obs` (cell metadata).
#
# We can start by identifying marker genes for each cluster. We can do this by setting `groupby='leiden'` to identify marker genes for each cluster defined by the `leiden` cluster assignments.
# %%
# via T-test
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#%% [markdown]
# `method=t-test` performs a Welch's t-test to identify genes that are differentially expressed between the groups. 
# 
# It does this by comparing the mean expression of each gene in each group in 1-vs-all comparisons. It then uses a Bonferroni correction to adjust the p-values for multiple comparisons.
#
# While the t-test is a relatively fast and simple test, it does make some assumptions about the data that may not be appropriate for single cell data. For example, it assumes that the data are normally distributed.
#
# We can also use `method='wilcoxon'` to perform a Wilcoxon rank sum test instead, which is a non-parametric test that is more appropriate for single cell data.
#%%
# Via Wilcoxon rank sum
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#%% [markdown]
# `rank_genes_groups` returns a structured array with the results of the test. 
# 
# We can access this using the `.uns` attribute of the anndata object.

#%%
adata.uns['rank_genes_groups']

#%% [markdown]
# Let's get a list of the ranked genes by cluster from this array.
# %%
ranked_genes_by_cluster = pd.DataFrame(adata.uns['rank_genes_groups']['names'])

ranked_genes_by_cluster.head(10)

#%% [markdown]
# This list can be exceptionally helpful in identifying cell types and states based on the expression of known marker genes.

#%% [markdown]
# ## Visualizing marker genes
# We can visualize the expression of these marker genes in a number of different ways.

# %%
sc.pl.rank_genes_groups_violin(adata, groups=["0"], n_genes=3 )

#%%
sc.pl.rank_genes_groups_heatmap(adata,n_genes=3)

#%%
sc.pl.rank_genes_groups_matrixplot(adata,n_genes=3,standard_scale='var')

# %%
sc.pl.violin(adata, ['Snap25','Slc17a6','Pax6','Gad1'], groupby='leiden')

#%% [markdown]
# We can get a list of the top 2 marker genes for each cluster and use this list to plot more information about this subset of genes.
# %%
marker_genes = ranked_genes_by_cluster.head(2).T.to_numpy().flatten().tolist()

marker_genes
#%%
sc.pl.dotplot(adata, marker_genes, groupby='leiden')

# %%
sc.pl.stacked_violin(adata, marker_genes, rotation=90, groupby='leiden')

#%% [markdown]
# Let's take this opportunity to write the results of this analysis to a file for future use.
#
# We will write the whole anndata `adata` object to a `.h5ad` file. 
# 
# This will allow us to easily load the results of this analysis into a new python session for further analysis or visualization.
#
# This can be done at any time with the `.write()` method of the anndata object.

# %%
adata.write("E18_mouse_cortex.h5ad")

# %% [markdown]
# Here's where we will end our introductory analysis.  
#
# The section below showcases an additional analysis using `scFates` to identify cell state transitions and pseudotime ordering of the subset of cells that are excitatory neurons.
#
# This is an optional section that is included for those who are interested in learning more about these topics.

# Subset to only excitatory neurons and their progenitors
#%% [markdown]
# ## Pseudotemporal analysis of Excitatory neuron development
#
# ### Subset Excitatory Neurons
# Identify clusters of excitatory neurons (main partition) and create a new anndata object with only these cells.

#%%
sc.pl.umap(adata,color='leiden',legend_loc='on data', add_outline=True)

ex_neuron_clusters = ["0","1","2","3","5","6","8","9"]

adata_ex = adata[adata.obs['leiden'].isin(ex_neuron_clusters),:].copy()

#%% [markdown]
# Import the scFates python module
#%%
import scFates as scf

#%% [markdown]
# Redo-PCA analysis on the subset of excitatory neurons

#%%
sc.tl.pca(adata_ex, svd_solver='arpack')

#%% 
sc.pl.pca(adata_ex,color=["Pax6"])

#%% [markdown]
# Fit a principal curve to the first two principal components

#%%
scf.tl.curve(adata_ex,Nodes=12,use_rep="X_pca",ndims_rep=2,)

#%%
scf.pl.graph(adata_ex,basis="pca")

#%% [markdown]
# Identify the root of the principal curve based on the expression of the radial glial marker gene 'Pax6' (marks progenitors).

# %%
scf.tl.root(adata_ex,"Pax6")

#%% [markdown]
# Calculate pseudotime position for each cell based on the principal curve
# %%
scf.tl.pseudotime(adata_ex,n_jobs=20,n_map=100,seed=42)

# %%
scf.pl.trajectory(adata_ex,basis="pca",arrows=True,arrow_offset=3)

#%% [markdown]
# Identify 'milestones' as key transitions along the trajectory
# %%
sc.pl.pca(adata_ex,color="milestones",legend_loc='on data', add_outline=True)

#%% [markdown]
# Rename milestones to provide more context
# %%
scf.tl.rename_milestones(adata_ex,new={"2":"Progenitors (RG)","0": "Maturing Neurons"})

# %%
scf.pl.milestones(adata_ex,basis="pca",annotate=True)

# %%
sc.pl.pca(adata_ex,color=["Pax6","Eomes","Snap25"])

#%% [markdown]
# Identify genes that are differentially expressed along the trajectory
#%%
scf.tl.test_association(adata_ex,n_jobs=20)

#%%
scf.tl.test_association(adata_ex,reapply_filters=True,A_cut=.5)

#%% [markdown]
# Fit a generalized additive model to the expression of each gene along the trajectory to make a smooth trend
#%%
scf.tl.fit(adata_ex,n_jobs=20)

#%% [markdown]
# Plot the expression of a few genes along the trajectory

# %%
scf.pl.single_trend(adata_ex,"Pax6",basis="pca",color_exp="k")
scf.pl.single_trend(adata_ex,"Eomes",basis="pca",color_exp="k")
scf.pl.single_trend(adata_ex,"Snap25",basis="pca",color_exp="k")

#%% [markdown]
# Cluster the genes based on their expression trends along the pseudotime trajectory

# %%
scf.tl.cluster(adata_ex,n_neighbors=100,metric="correlation")

#%% [markdown]
# Draw a heatmap of the expression of the genes in each cluster along the trajectory
#%%
for c in adata_ex.var["clusters"].unique():
    scf.pl.trends(adata_ex,features=adata_ex.var_names[adata_ex.var.clusters==c],basis="pca")


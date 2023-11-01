
#%%
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

#%% 
# Download dataset
matrix_url = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_filtered_feature_bc_matrix.tar.gz"
h5_url = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_filtered_feature_bc_matrix.h5"

#%%
# Some default settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%%
adata = sc.read_10x_h5('data/neuron_10k_v3_filtered_feature_bc_matrix.h5',gex_only=False)

# adata = sc.read_10x_mtx(
#     'data/filtered_feature_bc_matrix',  # the directory with the `.mtx` file
#     var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
#     cache=True)

#%%
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

#%%
adata

#%%
sc.pl.highest_expr_genes(adata, n_top=20, )

#%% [markdown]
### Basic filtering

#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#%%
# Identify mitochondrial genes by gene name
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'

#%%
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# %%
adata.var

#%% 
adata.obs

# %%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# %%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# %%
adata = adata[adata.obs.pct_counts_mt < 20, :]

# %%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# %%
# Total-count normalize (library size correct) the data matrix
sc.pp.normalize_total(adata, target_sum=1e4)

#%%
# Log transform
sc.pp.log1p(adata)

# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# %%
sc.pl.highly_variable_genes(adata)

# %%
adata.raw = adata

# %%
adata = adata[:, adata.var.highly_variable]

# %%
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#%%
sc.pp.scale(adata, max_value=10)

# %%

#%%
sc.tl.pca(adata, svd_solver='arpack')

#%%
sc.pl.pca(adata)


# %%
sc.pl.pca(adata, color="n_genes")
sc.pl.pca(adata, color="pct_counts_mt")

# %%
sc.pl.pca_variance_ratio(adata, log=True)

# %%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#%%
sc.tl.leiden(adata, resolution=0.3)

# %%
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph

#%%
sc.pl.paga(adata, color='leiden')

#%%
sc.pl.paga(adata, color=['Pax6','Eomes','Snap25'])

#%%
sc.tl.umap(adata, init_pos='paga')

# %%
sc.pl.umap(adata,color='leiden', add_outline=True)

#%%
sc.pl.umap(adata,color=['Pax6','Sox2','Top2a'], add_outline=True)

# %%
sc.pl.umap(adata,color=['Eomes','Tbr1','Slc17a6'], add_outline=True)

# %%
sc.pl.umap(adata,color=['Cux2','Rorb','Fezf2'], add_outline=True)

# %%
sc.pl.umap(adata,color=['Gad1','Gad2','Slc32a1'], add_outline=True)

# %%
# via T-test
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#%%
# Via Wilcoxon rank sum
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# %%
ranked_genes_by_cluster = pd.DataFrame(adata.uns['rank_genes_groups']['names'])

# %%
sc.pl.rank_genes_groups_violin(adata, groups=["0","1","2"], n_genes=3, )

# %%
sc.pl.violin(adata, ['Snap25','Slc17a6','Pax6','Gad1'], groupby='leiden')

# %%
marker_genes = ranked_genes_by_cluster.head(2).T.to_numpy().flatten().tolist()

sc.pl.dotplot(adata, marker_genes, groupby='leiden')

# %%
sc.pl.stacked_violin(adata, marker_genes, rotation=90, groupby='leiden')

# %%
adata.write("E18_mouse_cortex.h5ad")

# %%
# from sklearn.decomposition import NMF 

# model = NMF(n_components=20, init='random', random_state=0)

# W = model.fit_transform(adata.raw.X)
# H = model.components_

# Subset to only excitatory neurons and their progenitors
#%%
sc.pl.umap(adata,color='leiden')

ex_neuron_clusters = ["0","1","3","4","5","6","8","9"]

adata_ex = adata[adata.obs['leiden'].isin(ex_neuron_clusters),:].copy()



#%%
import pandas as pd

# %%
dat = pd.read_csv('data/rubstein_counts.txt',index_col=0)
# %%
dat.index

# %%
genes = pd.DataFrame(dat.index,columns=['gene_id'])
genes = genes.set_index('gene_id')
genes.to_csv('data/rubstein_genes.txt')

# %%

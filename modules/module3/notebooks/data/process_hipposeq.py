# Process Hippo-Seq data for AnnData input

#%%
import pandas as pd
import gget

#%%
dat = pd.read_csv('GSE74985_genes.read_group_tracking.txt',sep='\t')
dat_wide = dat.pivot_table(index='tracking_id',columns=['condition','replicate'],values='FPKM',aggfunc='mean')

#%%
gene_info = pd.DataFrame({'gene_id':dat_wide.index.to_list()})

# %%
# This is unbearably slow
#results = gget.info(gene_info['gene_id'].to_list(),verbose=True)
gencode_table = pd.read_csv('annotation.tbl',sep='\t',header=None)
gencode_table.columns = ['transcript_id', 'gene_id','chr', 'start', 'end', 'strand', 'exons', 'gene_biotype', 'gene_name']

gencode_subset = gencode_table[['gene_id','chr','start','end','strand','gene_biotype','gene_name']].copy()
gencode_subset_group = gencode_subset.groupby('gene_id')
#%%
genes = pd.DataFrame({'gene_id':gencode_subset['gene_id'].unique()})
genes['chr'] = gencode_subset_group['chr'].first().to_list()
genes['start'] = gencode_subset_group['start'].min().to_list()
genes['end'] = gencode_subset_group['end'].max().to_list()
genes['strand'] = gencode_subset_group['strand'].first().to_list()
genes['gene_biotype'] = gencode_subset_group['gene_biotype'].first().to_list()
genes['gene_name'] = gencode_subset_group['gene_name'].first().to_list()

genes['gene_id_short'] = genes['gene_id'].str.split('.',expand=True)[0]


#%%
# Left join gene_info and genes to get gene info
gene_info = gene_info.merge(genes,left_on='gene_id',right_on='gene_id_short',how='left')
gene_info = gene_info.rename(columns={"gene_id_x": "gene_id", "gene_id_y": "gene_id_version"})

#%%
assert(all(gene_info.gene_id == dat_wide.index.to_list()))


#%%
#gene_info[gene_info['gene_id'].str.startswith('ERCC')]['gene_biotype'] = 'ERCC_spikein'

#%%
gene_info.to_csv('GSE74985_gene_info.csv',index=False)

#############
# Sample info
#############
# %%
sample_info = pd.DataFrame({'sample':dat_wide.columns.map(lambda x: x[0]+'_'+str(x[1])).to_list()})

# %%
sample_info[['tissue','position']] = sample_info['sample'].str.split('_',expand=True)[[0,1]]
sample_info.loc[sample_info['position'].isin([str(x) for x in range(3)]), 'position'] = None
sample_info['replicate'] = [0,1,2]*8

# %%
sample_info.to_csv('GSE74985_sample_info.csv',index=False)

# %%
dat_wide.columns = sample_info['sample'].to_list()
dat_wide.index = gene_info['gene_id'].to_list()
dat_wide.to_csv('GSE74985_data.csv')
# %%

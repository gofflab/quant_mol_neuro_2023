#%%

# Install Tools

#%% [bash]
%%bash
mamba install -c bioconda hisat2 stringtie fastqc multiqc anndata

# Get the data

#%%
%%bash
mkdir data

cp /data/me440_lgoff2/datasets/RNA-Seq/data/raw/GSE74985/*.fastq.gz data/

#%% [markdown]
# # Fetch mouse reference genome GRCm29 (.fa) and reference annotation (.gtf)

#%%
%%bash
refDir="reference"
mkdir $refDir
genomeURL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz"
annotationURL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz"

wget --directory-prefix=$refDir $genomeURL
wget --directory-prefix=$refDir $annotationURL


#%% [bash]
%%bash
gunzip reference/*.gz

#%% [bash]

#%%
# Run fastqc

#%%
%%bash
mkdir fastqc_output

#%%
%%bash
for file in data/*.fastq.gz
do
	sbatch --job-name=$file-fastqc -o $file-fastqc.out -e $file-fastqc.err --wrap="fastqc -o fastqc_output $file"
done

#%%
# Build Hisat2 index
#%% 
%%bash
mkdir hisat2_index

#%%
%%bash
jobName="hisat_index"
sbatch --job-name=$jobName -o $jobName.out -e $jobName.err --ntasks=1 --cpus-per-task=16 --mem=12G --time=4:00:00 --wrap="hisat2-build -p 16 -f reference/GRCm39.primary_assembly.genome.fa hisat2_index/GRCm39.primary_assembly.genome"

# Run Hisat2 for each sample

#%% [bash]
%%bash
mkdir hisat2_alignments

#%% [bash]
%%bash
for file in data/*.fastq.gz
do
    base=$(basename $file .fastq.gz)
    sbatch --job-name=$base-hisat2 -o hisat2_alignments/$base-hisat2.out -e hisat2_alignments/$base-hisat2.err --cpus-per-task=16 --mem=12G --wrap="hisat2 -p 16 -x hisat2_index/GRCm39.primary_assembly.genome -U $file -S hisat2_alignments/$base.sam"
done

# Convert SAM to BAM

#%% [bash]
%%bash
for file in hisat2_alignments/*.sam
do
    base=$(basename $file .sam)
    sbatch --job-name=$base-sam2bam -o hisat2_alignments/$base-sam2bam.out -e hisat2_alignments/$base-sam2bam.err --wrap="samtools view -bS $file > hisat2_alignments/$base.bam"
done

# Sort BAM files by chromosome position
#%% [bash]
%%bash
for file in hisat2_alignments/*.bam
do
    base=$(basename $file .bam)
    sbatch --job-name=$base-sort -o hisat2_alignments/$base-sort.out -e hisat2_alignments/$base-sort.err --wrap="samtools sort $file -o hisat2_alignments/$base.sorted.bam"
done

# Clean up the SAM files
# Estimated time: <30 seconds
#%% [bash]
%%bash
rm hisat2_alignments/*.sam

# Index the BAM files
# Estimated time: 2-3 minutes
#%% [bash]
%%bash
for file in hisat2_alignments/*.sorted.bam
do
    base=$(basename $file .sorted.bam)
    sbatch --job-name=$base-index -o hisat2_alignments/$base-index.out -e hisat2_alignments/$base-index.err --wrap="samtools index $file"
done

# Run feature quantification with StringTie
#%% [bash]
%%bash
mkdir stringtie_quant

#%% [bash]
%%bash
for file in hisat2_alignments/*.sorted.bam
do
    base=$(basename $file .sorted.bam)
    sbatch --job-name=$base-stringtie -o stringtie_quant/$base-stringtie.out -e stringtie_quant/$base-stringtie.err --cpus-per-task=4 --wrap="stringtie -p 4 -e -G reference/gencode.vM33.annotation.gtf -A stringtie_quant/${base}_abundances.tab $file"
done

#%% [markdown]
### Run multiqc

#%% [bash]
%%bash
mkdir reports

#%% [bash]
%%bash
multiqc --outdir=reports .

# Let's aggregate all of the abundance files into a single pandas dataframe.  Now we're switching over to python.

#%%
import pandas as pd
from functools import reduce
import glob

#%%
# Make a list of all abundance files
abundance_files = glob.glob("stringtie_quant/*_abundances.tab")
abundance_files.sort() #In place sorting alphanumerically

# Parse abundance files to get sample names (in same order)
sample_names = [x.split("/")[1].split("_")[0] for x in abundance_files]

#%%
#group_value = ["Gene ID","Gene Name"]
group_value = ["Gene Name"]

# Use list comprehension to read in all abundance files (in order) into a list of pandas dataframes
abundances_list = [pd.read_csv(x, sep="\t",index_col=group_value) for x in abundance_files]

expr = pd.DataFrame({sample_names[0]:abundances_list[0]["TPM"]}, index=abundances_list[0].index)
expr = expr.groupby(group_value).sum("TPM") # Sum TPMs for fragmented genes output as duplicate rows by StringTie

for i in range(1,len(abundances_list)):
    expr_ = abundances_list[i]["TPM"]
    expr_ = expr_.groupby(group_value).sum("TPM")
    expr = expr.join(expr_, how="outer")
    expr.rename(columns={expr.columns[-1]: sample_names[i]},inplace=True)

#%%
# Check that the columns of the expression dataframe match sample names
all(expr.columns == sample_names)

#%%
# Grab the gene IDs from the index of the expression dataframe and populate a 'geneInfo' dataframe
geneInfo = pd.DataFrame(expr.index, columns=group_value)
geneInfo = geneInfo.set_index(group_value)


#%% [markdown]
#### Get the sample metadata
metadata_file = "GSE74985_sample_info.csv"

metadata = pd.read_csv(metadata_file, index_col=["Run"])

# Trim metadata to only those columns of interest
metadata = metadata[["cell_type","Location","Organism","Sample Name","source_name","tissue"]]

metadata[["position","region","cell"]] = metadata["source_name"].str.split(" ",2,expand=True)

metadata.head()

# %%
import anndata as ad

#%%
# Create AnnData object
adata = ad.AnnData(X=expr.T, obs=metadata, var=geneInfo)

# %%
# Remove genes with no expression in any sample
adata = adata[:,adata.X.sum(axis=0)>0]

#%%
import plotnine as pn

def plot_gene(adata ,gene_id):                       
    dat_ = adata[:,gene_id].copy()
    plot_df = dat_.obs.copy()
    plot_df["TPM"] = dat_.X.flatten()
    p =  (                                                      
        pn.ggplot(                                              
            plot_df,     
            pn.aes(x="source_name", y="TPM", fill="Location"),        
        )
        + pn.geom_boxplot(outlier_alpha=0.0)
        + pn.geom_point(size = 1)                               
        + pn.ggtitle(gene_id)                                   
        + pn.xlab("Sample")                          
        + pn.ylab("Gene Expression (TPM)")        
        + pn.labs(color="Location") 
        + pn.theme(legend_position="bottom")
    )
    return p  

# %%
plot_gene(adata,"Dll1")
# %%

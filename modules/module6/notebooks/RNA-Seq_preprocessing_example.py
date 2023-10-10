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
    sbatch --job-name=$base-stringtie -o stringtie_quant/$base-stringtie.out -e stringtie_quant/$base-stringtie.err --cpus-per-task=16 --wrap="stringtie -p 16 -e -G reference/gencode.vM33.annotation.gtf -A stringtie_quant/${base}_abundances.tab $file"
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
# Use list comprehension to read in all abundance files (in order) into a list of pandas dataframes
abundances_list = [pd.read_csv(x, sep="\t",index_col=0) for x in abundance_files]

for x in abundances_list:
    x.index = x.index.str.strip()

abundances_list = [x.sort_index() for x in abundances_list]


#%%
# Grab the gene information columns from the first dataframe
#geneInfo = pd.DataFrame(abundances_list[0][["Gene ID","Gene Name","Reference","Strand","Start","End"]])

# Concatenate the TPM columns from all dataframes into a single dataframe of expression estimates (tpm)
#expr = pd.concat([x["TPM"] for x in abundances_list], verify_integrity=True, axis=1)

expr = pd.DataFrame({sample_names[0]:abundances_list[0]["TPM"]}, index=abundances_list[0].index)

for i in range(1,len(abundances_list)):
    expr = expr.join(abundances_list[i]["TPM"], how="outer")
    expr.rename(columns={expr.columns[-1]: sample_names[i]},inplace=True)

expr.shape


# Rename the columns of the expression dataframe to the sample names
expr.columns = sample_names

# Concatenate the gene information and expression dataframes into a single dataframe for viewing
dat = pd.concat([geneInfo,expr], axis=1)

#%%
dat.head()

#%% [markdown]
#### Get the sample metadata
metadata_file = "GSE74985_sample_info.csv"

metadata = pd.read_csv(metadata_file, index_col=0)

# Trim metadata to only those columns of interest
metadata = metadata[["cell_type","Location","Organism","Sample Name","source_name","tissue"]]

metadata.head()

# %%
import anndata as ad

#%%
# Create AnnData object
adata = ad.AnnData(X=expr.T, obs=metadata, var=geneInfo)
# %%

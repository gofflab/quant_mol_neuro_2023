# %% [markdown]
# # Module 3 - Sequencing Data Formats and Structures & Intro to Classes and Objects

# %% [markdown]
# ## Sequencing Data Formats and Structures
# Modern sequencing technologies generate large amounts of data. These data is usually stored in files that are formatted in a specific way. 
# Most of these file formats are designed to be easily parsed by computers, _and_ to be easily read by humans. 
#
# There are different file formats for different types of genomic data. For example, there are different file formats for storing nucleotide sequences, the mapping (or 'alignment') of sequences to a reference, and for annotations such as genomic features (e.g. genes, transcripts, exons, etc).
#
# Each file format has a specific structure, and contains specific information.
#
# In this section, we will look at some of the more common file formats used in genomics, and how to read and parse them using python. 

#%% [markdown]
# ### File formats for storing sequence information
# There are several different file formats for storing nucleotide or polypeptide sequence information. 
#
# #### FASTA
# The [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file format is a simple, text-based format for representing either nucleotide sequences or peptide sequences, in which nucleotides or amino acids are represented using [single-letter codes](https://www.bioinformatics.org/sms/iupac.html).
# 
# A sequence in FASTA format begins with a single-line description, followed by lines of sequence data. The description line is distinguished from the sequence data by a greater-than (`>`) symbol at the beginning. 
# 
# An example sequence in FASTA format is:
#
# ```
# >sequence1
# ATGCAGTACTGACGTATCGCATTCGTCATGC
# ```
#
# The first line in this example is the description, and the following line is the sequence itself. The description line is always the first line of a record, and must begin with a `>` symbol in the first column. 
# 
# The string following the `>` symbol is the identifier (name) of the sequence, and the rest of the line is an optional description of the entry. 
# 
# There should be no space between the `>` and the first letter of the identifier. 
# 
# The sequence ends if another line starting with a `>` appears; this indicates the start of another sequence.
#
# A sample .fasta file can be found in the `data` directory of this module as `data/sample_sequences.fa`.

#%% [markdown]
# ### FASTQ
# The [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file format is similar to the `FASTA` format, but it also contains [quality scores](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm) for each base derived from the sequencing instrument on which a given sequence was generated.
#
# This is a common format for 'raw' sequencing data, and the per-base scores are used downstream to help inform the alignment of the reads to a reference sequence (e.g. genome/transcriptome).
#
# The quality scores are usually represented as ASCII characters, with each character representing a calculated numeric score. 
# 
# The numeric scores are calculated by the sequencing instrument, and are usually a measure of the probability that the base call is incorrect.
# 
# `FASTQ` formatted sequences also begin with a single-line description, followed by lines of sequence data. The description line is distinguished from the sequence data by an `@` symbol at the beginning. Followed by the sequence.
#
# After the sequence, there is another description line, which is distinguished from the quality string by a plus (`+`) symbol at the beginning.
#
# Finally, the per-base quality string is included on the next line. This line must be the same length as the sequence itself, and the first character of the quality string corresponds to the first base of the sequence.
#
# An example sequence in `FASTQ` format is:
# ```
# @SEQ_ID
# GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
# +
# !''*((((***+)) % % % ++)( % %%%).1***-+*''))**55CCF>>>>>>CCC
# ```


#%% [markdown]
## File formats for storing read alignment information

#%% [markdown]
# ### SAM/BAM
# The [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) (Sequence Alignment/Map) format is a generic format for storing large nucleotide sequence _**alignments**_. SAM is a text format, and is human readable. 
# 
# SAM files, and their compressed binary equivalent, [BAM](https://en.wikipedia.org/wiki/SAM_(file_format)#Binary_Alignment_Map_(BAM)), are currently the standard output format for most sequence alignment programs.
#
# SAM files are tab-delimited, with each line in the file representing a single read alignment to a reference.
# 
# A SAM file has two main parts, a header section and an alignment section. 
# 
# The header section contains information about the reference sequence(s) that the reads were aligned to. Header rows start with the `@` symbol and provide details about the reference sequences, the alignment program used, and other information about the alignment.
# 
# An example SAM file header:
# ```
# @HD	VN:1.6	SO:unsorted
# @SQ	SN:chr1	LN:248956422
# @SQ	SN:chr2	LN:242193529
# @SQ	SN:chr3	LN:198295559
# @PG	ID:bwa	PN:bwa	VN:0.7.12
# ```
#
# Some details about the SAM format from the above example:
# - @HD: Header line indicating format version.
# - @SQ: Sequence dictionary for reference sequences. In this case, we have reference sequences for chr1 and chr2.
# - @PG: Program record. Here, it indicates that the BWA tool, version 0.7.12 was used to create the alignment.
#
# For a detailed description of the SAM file header, see the [SAM file format specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

#%% [markdown]
# The alignment section begins immediately after the header and contains the actual read alignment information. 
# 
# Each line in the alignment section represents a single read alignment to a reference sequence.
#
# Example SAM file record:
# ```
# SRR001666.1	163	chr1	10000	255	4M	*	0	0	*	*
# ```
# 
# The first 11 columns of a SAM/BAM/CRAM file record are the same, and contain the following information:
#
# | Column | Description | Details |
# | --- | --- | --- |
# | 1 | QNAME - Query template NAME | e.g. read name|
# | 2 | FLAG - bitwise FLAG | a bitflag that encodes information about the alignment |
# | 3 | RNAME - Reference sequence NAME | e.g. chromosome name |
# | 4 | POS - 1-based leftmost mapping POSition | the 'start' position of the alignment |
# | 5 | MAPQ - MAPping Quality | a phred-scaled probability that the alignment is incorrect |
# | 6 | CIGAR - CIGAR string | a string encoding the specific details of the alignment |
# | 7 | RNEXT - Ref. name of the mate/next read | used for paired-end alignments |
# | 8 | PNEXT - Position of the mate/next read | used for paired-end alignments |
# | 9 | TLEN - observed Template LENgth | a measure of the observed length of the fragment (insert size) |
# | 10 | SEQ - segment SEQuence | the actual sequence of the aligned portion of the read |
# | 11 | QUAL - ASCII of Phred-scaled base QUALity+33 ||
#
# The remaining columns are optional, and contain additional information about the alignment. These optional columns are often referred to as 'tags', and are used to store additional information about the alignment.
#
# The [samtools](http://www.htslib.org/) command-line tool is a popular tool for working with SAM/BAM files. It can be used to convert between SAM and BAM formats, sort and index BAM files, and extract specific information from BAM files.

#%% [markdown]
### CRAM
# An emerging standard for storing sequence alignments is the [CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format)) (Compressed Read Alignment Map) format. CRAM is a compressed binary format that is designed to be a more efficient alternative to the BAM format.
# 
# CRAM files are not as common as BAM files, but are becoming more popular as the standard for storing sequence alignments.
#
# CRAM files are also tab-delimited, with each line in the file representing a single read alignment to a reference.

### Other common alignment formats
# There are several other file formats for storing sequence alignments, including:
# - [MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) (Multiple Alignment Format)
# - [PSL](https://genome.ucsc.edu/FAQ/FAQformat.html#format2) (PSL format)
# - [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (Variant Call Format)
# - [WIG](https://genome.ucsc.edu/goldenPath/help/wiggle.html) (Wiggle format)
# - [BigWig](https://genome.ucsc.edu/FAQ/FAQformat.html#format9) (Big Wig format)
# - [BigBed](https://genome.ucsc.edu/FAQ/FAQformat.html#format8) (Big Bed format)

#%% [markdown]
# ## File formats for storing genomic feature information
#
# Another common type of file format in genomics is the file format for storing genomic feature information such as genes, transcripts, exons, promoters, alignment peaks, etc. 
# 
# These files describe the genomic features in a genome, and are often used to annotate sequence alignments or other genomic data.
#
# An important caveat for these file formats is that often, most of the information presented is 'relative' to a specific reference genome or sequence. For example, the chromosome, start, and end positions of a gene are specific to a particular reference genome and assembly.
#

#%% [markdown]
# ### BED
# The [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file format is a text-based format for representing genomic features (e.g. genes, transcripts, exons, peaks, or other 'intervals' along a genome).
# 
# BED files are tab-delimited text files, with each line in the file representing a single genomic feature. 
#
# BED files have three required fields, and six additional optional fields. The required fields are:
# - Chromosome
# - Start position
# - End position
#
# The optional fields are:
# - Feature Name
# - Score (user defined)
# - Strand
# - Thick start (for features that are not uniform width. e.g. exons)
# - Thick end 
# - Item RGB (used to colorize items)
#
# The following is an example BED file record for a gene (DDX11L1) on chromosome 1:
# ```
# chr1	11873	14409	DDX11L1	0	+	11873	11873	0	3	354,109,1189	0,739,3479
# ```
#
# BED files are often used to draw features on a genome browser, such as the [UCSC Genome Browser](https://genome.ucsc.edu/). Most of the tracks that encode annotated features on the genome browser are in BED format.

#%% [markdown]
# ### GFF/GTF
# The [GFF](https://genome.ucsc.edu/FAQ/FAQformat.html#format3) (General Feature Format) and the more specialized [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) (Gene Transfer Format) file formats are text-based formats for representing more complex genomic features (e.g. genes, transcripts, exons, etc).
# 
# These formats are similar to the BED format, but are more complex and allow for more information to be stored about each feature. 
# 
# For example, the GFF/GTF format allows for hierarchical relationships between features (e.g. exons are part of specific transcript isoforms, which belong to specific genes).
#
# Most transcriptome annotation files are provided in GFF/GTF format.
#
# GTF records have 9 required fields:
# - Chromosome (or other reference sequence identifier)
# - Source (e.g. the name of the program that generated this feature)
# - Feature type (usually one of gene, transcript, CDS, 5UTR, 3UTR, start_codon, stop_codon, or exon)
# - Start position
# - End position
# - Score (if applicable, otherwise '.')
# - Strand
# - Frame (if applicable, otherwise '.')
# - Attributes (a semicolon-separated list of tag-value pairs providing additional information about each feature)
#
# The Attributes field is specific to GTF files, and is used to store additional information about each feature
#
# The following is an example GTF file record for a gene (DDX11L1) on chromosome 1:
# ```
# chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
# ```
# 
# GTF files can include information about isoforms, exons, and other features. 
# 
# The following is an example GTF file record for an isoform of the DDX11L1 gene and it's 3 exons:
# ```
# chr1  HAVANA transcript   11869   14409   .   +   .   gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
# chr1  HAVANA  exon    11869   12227   .   +   .   gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
# chr1  HAVANA  exon    12613   12721   .   +   .   gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
# chr1  HAVANA  exon    13221   14409   .   +   .   gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 3; exon_id "ENSE00002312635.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
# ```
#
# 

#%% [markdown]
# The Gencode project ([gencodegenes.org](https://www.gencodegenes.org/)) is the current standard for gene annotation.
#
# They provide detailed, comprehensive gene annotation files in GTF format for the [human](https://www.gencodegenes.org/human/) and [mouse](https://www.gencodegenes.org/mouse/) genomes.
# 
# The Gencode annotation files contain information about genes, transcripts, exons, and other genomic features.
#
# Let's download the Gencode annotation file for the most recent release (v44) which uses the GRCh38.p14 human genome assembly.

#%% 
gencode_gtf = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz'

# Download the file and save it to the data directory
import urllib.request
import gzip
urllib.request.urlretrieve(gencode_gtf, 'data/gencode.v44.annotation.gtf.gz')

# %% [markdown]
# Uncompress the file `data/gencode.v44.annotation.gtf.gz` and open it in VSCode to take a look at the contents and familiarize yourself with the format.

#%% [markdown]
# ### Parsing and handling genomic data file formats
# Often times, we use tools that will work directly on these file formats. For example, we might use a tool to align sequencing reads to a reference genome, and the inputs and outputs of that tool will be in one of these file formats.
# 
# There are, however, times when we will need to examine or iterate through these files ourselves. For example, we might need to extract specific information from a file, or we might need to combine information from multiple files. 
#
# In order to do this, it's useful to know how to read and parse these files using python.
#
# Several Python libraries and tools have been developed to parse and handle genomic data file formats. Here's a rundown of some of the popular ones:
#
# * FASTQ:
# 	- [Biopython](https://biopython.org/)'s [SeqIO](https://biopython.org/wiki/SeqIO): The `SeqIO` module of Biopython provides a simple interface for reading from and writing to various bioinformatics file formats, including FASTQ.
# * FASTA:
# 	- [Biopython](https://biopython.org/)'s [SeqIO](https://biopython.org/wiki/SeqIO): As with FASTQ, the SeqIO module is versatile and can be used for reading and writing FASTA files.
# 	- [Pyfastx](https://pypi.org/project/pyfastx/): A a python module for fast parsing and access to FASTA/Q files.
# * SAM/BAM:
# 	- [Pysam](https://pysam.readthedocs.io/en/latest/): A Python module for reading and manipulating Samfiles. It's a wrapper around the samtools suite, which provides functions to read, manipulate, and write SAM/BAM format files.
# 	- [Bamnostic](https://bamnostic.readthedocs.io/en/latest/): It is an OS-agnostic library for working with BAM(Binary Alignment Map) files.
# * BED:
# 	- [pyBedTools](https://daler.github.io/pybedtools/): This is a wrapper around the popular BedTools suite, allowing for manipulation of genomic intervals and any associated metadata in a variety of formats(including BED).
# 	- Pandas: While not specific to genomics, Pandas can be used to handle BED files as simple tables, especially when there's a need for data manipulation tasks.
# * GFF/GTF:
# 	- [gffutils](https://daler.github.io/gffutils/): This is a tool for working with GFF and GTF format files, offering database creation, querying, and other operations on these formats.
# 	- Bcbio-gff: A Python library for reading and writing GFF(version 2 and 3) files.
#	- [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread): gffread is not a python library but actually a compiled program that can be used to convert GFF/GTF files to other formats, including BED, GTF, and FASTA. It is a fast and powerful tool for manipulating GFF/GTF files.
# 

# %% [markdown]
# ## Classes and objects
# Python is an object-oriented programming language. This means that it provides features that support object-oriented programming (OOP). 
# 
# Object-oriented programming brings together data and its behaviour (methods) in a single location (called an “object”) making it easier to conceptualize and understand. 
# 
# This reduces complexity and makes it easier to reuse code in different parts of a program or in different programs.
# 
# A `class` is a blueprint for an `object`. It describes how the `object` is made (what data it contains and what methods it has). An `object` is an instance of a `class`. It contains real values stored in locations defined by the `class`. 
# 
# You can create as many objects as you want from a `class`. Each `object` is independent of the others. You can modify an `object` without affecting the others.
# 
# ### Creating a class
# Lets start by creating a class to hold a DNA sequence. We will need to define the `attributes` of the class, which are the variables that will be associated with each instance of the class.
# 
# We will start with a simple class that only has a single attribute, the DNA sequence itself.

# %%
class DNA:
    def __init__(self, seq, name):
        self.seq = seq
        self.name = name

# %% [markdown]
# Now let's create an instance of our class. We do this by calling the class name as if it were a function, and passing the required arguments to the class `__init__` method. 
# 
# The `__init__` method is a special method that is called when an instance of the class is created. It is used to initialize the `attributes` of the class. 
# 
# The first argument of the `__init__` method (or _**any**_ method actually) is always a reference to the object itself (the instance of the class). By convention, this argument is called `self`. 
# 
# The other arguments are the ones that we passed to the class when we created the instance.

# %%
myDNA = DNA('ATGCAGTACTGACGTATCGCATTCGTCATGC', 'Sequence1')

# %% [markdown]
# Right now this instance doesn't do much. It has two attributes (`.seq` and `.name`), but no methods to do anything with it. 
# 
# Let's add a few methods to our class that will allow us to calculate features of the DNA sequence.
#
# Methods are created by defining functions within the class. As with the `__init__` methood, the first argument of any method is always a reference to the object itself (`self`).

# %%
class DNA:
    def __init__(self, seq, name):
        self.seq = seq
        self.name = name
    
    def __len__(self):
        return(len(self.seq))
    
    def gc(self):
        G = self.seq.count("G")
        C = self.seq.count("C")
        nGC = G+C
        return(nGC/len(self.seq)*100)
    
    def tm(self):
        G = self.seq.count("G")
        C = self.seq.count("C")
        A = self.seq.count("A")
        T = self.seq.count("T")
        Tm = 64.9+41*(G+C-16.4)/(A+T+G+C)
        return(Tm)
    
    def toFASTA(self):
        return(f">{self.name}\n{self.seq}")

#%% [markdown]
# We've seen some of these functions before in previous modules. Previously, we've used functions like this to operate on strings. Now, we've defined these functions as methods of our class, and they are now `owned` by all instances of the `DNA` class.
#
# Now we can create an instance of our class and use the methods to calculate specific features of the DNA sequence contained within the object.
# %%
myDNA = DNA('ATGCAGTACTGACGTATCGCATTCGTCATGC', 'Sequence1')

print(myDNA.gc())
print(myDNA.tm())

#%% [markdown]
# We also added a method to access the special python `len()` function for our class. This allows us to use `len()` on our class instances to return the length of a _**specific attribute**_ (in this case, the sequence length), and have it return the length of that attribute.

#%%
print(len(myDNA))

#%% [markdown]
# Finally, we added a method to return the sequence in FASTA format. This method returns a string that is formatted as a FASTA sequence. We can use this method to quickly 'dump' (print) the sequence in the object in the specific FASTA format.

#%% 
print(myDNA.toFASTA())

# %% [markdown]
# 

# %% [markdown]
# # Reading and parsing Genomics data files
# ## Reading FASTA files
# As outlined above, FASTA is a file format for representing nucleotide or peptide sequences. A FASTA file consists of a header line followed by lines of sequence data. The header line is distinguished from the sequence data by a greater-than (">") symbol in the first column. The word following the ">" symbol is the identifier (name) of the sequence, and the rest of the line is an optional description of the entry. There should be no space between the ">" and the first letter of the identifier. The sequence ends if another line starting with a ">" appears; this indicates the start of another sequence.
# 
# FASTA is a common format in bioinformatics for storing sequence strings. It is a simple format that is easy to parse. 
# 
# Let's write a function that reads a FASTA file and returns a list of `DNA` objects.

# %%
# Create a python function to take a fasta filename as an argument and return a list of DNA objects
def parse_fasta(filename):
    # Open the file
    fastaFile = open(filename, 'r')
    # Create an empty list to store the instances of the DNA sequence class
    sequences = []
    # Create an empty string to store the current sequence
    currentSequence = ''
    # Create an empty string to store the current sequence name
    currentName = ''
    # Loop through the lines in the file
    for line in fastaFile:
        # If the line starts with a >, we have a new sequence
        if line.startswith('>'):
            # If we have a current sequence, add it to the list
            if currentSequence != '':
                sequences.append(DNA(currentSequence, currentName))
            # reset the current sequence
            currentSequence = ''
            # Set the current name to the name on the current header line
            currentName = line.strip()[1:]
        else:
            # Add the line to the current sequence
            currentSequence += line.strip()
    # Add the last sequence to the list
    sequences.append(DNA(currentSequence, currentName))
    # Close the file
    fastaFile.close()
    # Return the list of sequences
    return sequences

# %% [markdown]
# Now we can use our function to read a FASTA file and return a list of sequences.
#
# Let's test it out on the FASTA file containing some sample DNA sequences.

#%%
fasta_file = 'data/sample_sequences.fa'

sequences = parse_fasta(fasta_file)

sequences

#%%
# Extract the sequences from the list using the `.seq` atrribute 
[x.seq for x in sequences]

#%%
[x.gc() for x in sequences]

#%%
[print(x.toFASTA()) for x in sequences]
# %% [markdown]
# ## AnnData class and objects
# One of the more useful data formats in python for gene expression data is the [`AnnData`](https://anndata.readthedocs.io/en/latest/) object. 
#
# While originally developed to store and manipulate single-cell RNA-seq data, the `AnnData` class is a very flexible data structure that can be used to store and manipulate any type of gene expression or other genomic data.
# 
# `AnnData` is a `class` that is part of the `anndata` package and is used to store and manipulate gene expression data.
#
# At it's core, an `AnnData` object is a matrix of gene expression values (`X`). It also contains information about the features (genes) and samples that are in the matrix.
#
# The [`AnnData` object](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html) is a complex class, with many attributes and methods. We will not go into all of the details here, but we will look at some of the more common attributes and methods that are used to access and manipulate the data in the object.
#
# Let's start by importing the `anndata` package and creating a simple `AnnData` object from a gene expression matrix.
#
# To create an `anndata` object, we need a few minimal pieces of information:
# - The gene expression matrix `X` (a `numpy` ndarray)
# - The sample information `obs` (a `list`, `dict`, or `pandas` `Series` or `DataFrame`)
# - The gene information `var` (a `list`, `dict`, or `pandas` `Series` or `DataFrame`)

#%%
import anndata as ad
import numpy as np

adata = ad.AnnData(X=np.array([[1, 2], [3, 4], [5, 6], [7, 8]]), 
                   obs={'obs_names': ['sample1', 'sample2','sample3','sample4'], 'condition': ['wildtype', 'wildtype','knockout','knockout']}, 
                   var={'var_names': ['gene1', 'gene2'], 'gene_type': ['protein_coding', 'protein_coding'], 'chromosome': ['chr1', 'chr2']})

adata
# %% [markdown]
# Now that we have an `AnnData` object instantiated, let's look at some of the attributes and methods that we can use to access and manipulate the data in the object.
#
# We can access the gene expression matrix `X` using the `.X` attribute of the object.
# %% 
adata.X

#%% [markdown]
# We can access the sample information `obs` using the `.obs` attribute of the object. This will return a `pandas` `DataFrame` containing the sample information in the same order as the _**rows**_ of the gene expression matrix `X`.
#
#%% 
adata.obs

#%% [markdown]
# We can access the gene information `var` using the `.var` attribute of the object. This will return a `pandas` `DataFrame` containing the gene information in the same order as the _**columns**_ of the gene expression matrix `X`.
#%%
adata.var

# %% [markdown]
# `AnnData` objects can be subsetted using the `[]` notation that we have previously used to subset `pandas` `DataFrame` objects.

#%%
# Subset the object to only include the first two samples (rows)
adata[:2,:]

#%%
# Subset the object to only include the first two genes (columns)
adata[:,:2]

# %%
# Subset using gene or sample names
adata[['sample2'], ['gene2']]

# %%
# Subset using attribute values
adata[adata.obs['condition']=='wildtype', adata.var['gene_type']=='protein_coding']

#%% [markdown]
# From each of these subset 'views' of the original object, we still have acess to all of the attributes and methods of the original object.
#
# Importantly, `AnnData` objects are designed so that when you subset or filter in this way, all of the accessory information (sample and gene information) is also subsetted or filtered to match the subsetted gene expression matrix `X`.

#%%
adata[adata.obs['condition']=='wildtype',:].obs

#%%
adata[adata.obs['condition']=='wildtype',:].var

#%%
adata[adata.obs['condition']=='wildtype',:].X

#%% [markdown]
# Let's create an `AnnData` object from a real-world gene expression dataset.
#
# We will use the Hippo-Seq bulk RNA-Seq dataset from last week's module.
#
# - [GSE74985](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74985)
# - Cembrowski MS, Wang L, Sugino K, Shields BC et al. Hipposeq: a comprehensive RNA-seq database of gene expression in hippocampal principal neurons. Elife 2016 Apr 26;5:e14997. PMID: [27113915](https://www.ncbi.nlm.nih.gov/pubmed/27113915)
#
# Let's start by reading in the three files that we need to create the `AnnData` object:
# - The gene expression matrix `X` (a `pandas` `DataFrame`)
# - The sample information `obs` (a `pandas` `DataFrame`)
# - The gene information `var` (a `pandas` `DataFrame`)

#%% 
import pandas as pd
expression = pd.read_csv('data/GSE74985_data.csv', index_col=0, header=0)
expression = expression.T # We need to transpose the 'genes X samples' matrix to a 'samples X genes' matrix for AnnData

# %%
gene_info = pd.read_csv('data/GSE74985_gene_info.csv', header=0)
gene_info.index = gene_info['gene_id'].to_list()
gene_info

# %%
sample_info = pd.read_csv('data/GSE74985_sample_info.csv', header=0)
sample_info.index = sample_info['sample'].to_list()
sample_info

# %% [markdown]
# Now we can create our `AnnData` object with these three matrices
#

#%%
adata = ad.AnnData(X=expression, 
                   obs=sample_info,
                   var=gene_info)

# %% [markdown]
# Now we can use the same methods and attributes that we used above to access and manipulate the data in the object.

#%%
adata.X # The gene expression matrix (now stored as a numpy ndarray)

#%%
adata.obs # The sample information (stored as a pandas DataFrame)

#%%
adata.var # The gene information (stored as a pandas DataFrame)

#%% [markdown]
# And we can use the attributes in the `.var` and `.obs` dataframes to subset the object as well.

#%%
adata[adata.obs['tissue']=='ca1',:].X # Subset the object to only include samples from the CA1 region and return the expression matrix (X).

# %%
adata[:,adata.var['gene_biotype'] == 'lncRNA']

# %% [markdown]
# A final note about subsetting `AnnData` objects. The default when subsetting is to return a _**view**_ of the original object. 
# 
# This means that the subsetted object is not a copy of the original object, but rather a 'pointer' to the original object.
# 
# Any changes that are made to a view of the object will also be made to the original object itself.
#
# To avoid this, and to create a copy of the subsetted data, you can use the `.copy()` method of the object.

#%%
adata_lncRNA = adata[:,adata.var['gene_biotype'] == 'lncRNA'].copy()

adata_lncRNA
# %% [markdown]
# # Module 2 - Python classes, modules, and packages

# %% [markdown]
# ## Sequencing Data Formats and Structures

#%% [markdown]
# ### FASTA
# The [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file format is a text-based format for representing either nucleotide sequences or peptide sequences, in which nucleotides or amino acids are represented using single-letter codes. 
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
# The first line in this example is the description, and the following line is the sequence itself. The description line is always the first line of a record, and must begin with a greater-than symbol in the first column. 
# 
# The string following the `>` symbol is the identifier (name) of the sequence, and the rest of the line is an optional description of the entry. 
# 
# There should be no space between the `>` and the first letter of the identifier. 
# 
# The sequence ends if another line starting with a `>` appears; this indicates the start of another sequence.

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
# After the sequence, there is another description line, which is distinguished from the quality string by a plus (`+`) symbol at the beginning. This line is optional, and is usually used to indicate the sequencing instrument used to generate the sequence.
#
# Finally, the per-base quality string is included on the next line. This line must be the same length as the sequence itself, and the first character of the quality string corresponds to the first base of the sequence.
#
# An example sequence in `FASTQ` format is:
# ```
# @sequence1
# ATGCAGTACTGACGTATCGCATTCGTCATGC
# +
# ?A>EFF;FGGGGGGGGGGGGGGGGGGGGGGG
# ```


#%% [markdown]
# ### SAM/BAM
# The [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) (Sequence Alignment/Map) format is a generic format for storing large nucleotide sequence _**alignments**_. SAM is a text format, and is human readable. 
# 
# SAM files, and their compressed binary equivalent, [BAM](https://en.wikipedia.org/wiki/SAM_(file_format)#Binary_Alignment_Map_(BAM)), are currently the standard output format for most sequence alignment programs.
#
# SAM files are tab-delimited, with each line in the file representing a single read alignment to a reference.
# 
# Example SAM file record:
# ```
# SRR001666.1	163	chr1	10000	255	4M	*	0	0	*	*
# ```

#%% [markdown]
# ### BED
# The BED file format is a text-based format for representing genomic features (e.g. genes, transcripts, exons, peaks, or other 'intervals' along a genome).
# 
# 

#%% [markdown]
# ### GFF/GTF
# The GFF (General Feature Format) and GTF (Gene Transfer Format) file formats are text-based formats for representing more complex genomic features (e.g. genes, transcripts, exons, etc). 
# 
# These formats are similar to the BED format, but are more complex and allow for more information to be stored about each feature. For example, the GFF/GTF format allows for hierarchical relationships between features (e.g. exons are part of specific transcript isoforms, which belong to specific genes).
#
# 


#%% [markdown]
#

# %% [markdown]
# ## Classes and objects
# Python is an object-oriented programming language. This means that it provides features that support object-oriented programming (OOP). Object-oriented programming brings together data and its behaviour (methods) in a single location (called an “object”) making it easier to conceptualize and understand. 
# 
# This reduces complexity and makes it easier to reuse code in different parts of a program or in different programs.
# 
# A `class` is a blueprint for an `object`. It describes how the `object` is made (what data it contains and what methods it has). An `object` is an instance of a `class`. It contains real values stored in locations defined by the `class`. 
# 
# You can create as many objects as you want from a `class`. Each `object` is independent of the others. You can modify an `object` without affecting the others.
# 
# ### Creating a class
# Lets start by creating a class to hold a DNA sequence. We will need to define the attributes of the class, which are the variables that will be associated with each instance of the class.
# 
# We will start with a simple class that only has a single attribute, the DNA sequence itself.

# %%
class DNA:
	def __init__(self, seq, name):
		self.seq = seq
		self.name = name

# %% [markdown]
# Now let's create an instance of our class. We do this by calling the class name as if it were a function, and passing the required arguments to the class `__init__` method. The `__init__` method is a special method that is called when an instance of the class is created. It is used to initialize the attributes of the class. The first argument of the `__init__` method is always the object itself (the instance of the class). By convention, this argument is called `self`. The other arguments are the ones that we passed to the class when we created the instance.

# %%
myDNA = DNA('ATGCAGTACTGACGTATCGCATTCGTCATGC', 'Sequence1')

# %% [markdown]
# Right now this instance doesn't do much. It has two attributes (`.seq` and `.name`), but no methods to do anything with it. Let's add a few methods to our class that will allow us to calculate features of the DNA sequence.

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
# FASTA is a file format for representing nucleotide or peptide sequences. A FASTA file consists of a header line followed by lines of sequence data. The header line is distinguished from the sequence data by a greater-than (">") symbol in the first column. The word following the ">" symbol is the identifier (name) of the sequence, and the rest of the line is an optional description of the entry. There should be no space between the ">" and the first letter of the identifier. The sequence ends if another line starting with a ">" appears; this indicates the start of another sequence.
# 
# FASTA is a common format in bioinformatics for storing sequence strings. It is a simple format that is easy to parse. Let's write a function that reads a FASTA file and returns a list of DNA sequences.

# %%
# Create a python function to take a fasta filename as an argument an return a list of sequences
def parse_fasta(filename):
	# Open the file
	fastaFile = open(filename, 'r')
	# Create an empty list to store the sequences
	sequences = []
	# Create an empty string to store the current sequence
	currentSequence = ''
	# Loop through the lines in the file
	for line in fastaFile:
		# If the line starts with a >, we have a new sequence
		if line.startswith('>'):
			# If we have a current sequence, add it to the list
			if currentSequence != '':
				sequences.append(currentSequence)
			# Reset the current sequence
			currentSequence = ''
		# Otherwise, we have a sequence line
		else:
			# Add the line to the current sequence
			currentSequence += line.strip()
	# Add the last sequence to the list
	sequences.append(currentSequence)
	# Close the file
	fastaFile.close()
	# Return the list of sequences
	return sequences

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
			# Reset the current sequence and name
   			currentSequence = ''
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
# Let's test it out on a FASTA file containing some DNA sequences.

	   
	


# %% [markdown]
# ## AnnData class and objects
# One of the more common data formats in python for gene expression data is the [`AnnData`](https://anndata.readthedocs.io/en/latest/) object. 
# 
# This is a class that is part of the `anndata` package and is used to store and manipulate gene expression data.
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

adata = ad.AnnData(X=np.array([[1, 2], [3, 4]]), 
                   obs={'obs_names': ['cell1', 'cell2']}, 
                   var={'var_names': ['gene1', 'gene2']})

adata

# %% 
adata.X

#%% 
adata.obs

#%%
adata.var

#%% [markdown]
# Let's create an `AnnData` object from a gene expression matrix.

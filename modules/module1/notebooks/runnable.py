"""
Similar to how we provided arguments to functions, we can also flexibly provide python scripts with arguments.

This is an example code.

The syntax in bash for this is:
$ python runnable.py ACTGATGCTAGCTGACTGATCTAGCTGA
> TDAS_LI_L

$ python runnable.py ACTGATGCTAGCTGACTGATCTAGCTGA --ignore_stop
> TDASLIL
"""

import click

# fmt: off
GENCODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}
# fmt: on


# First let's write a function to return the amino acid associated with a given codon.
def lookup_AA(codon: str) -> str:
    # get the value associated with the codon key, if not found, return "X"
    aa = GENCODE.get(codon, "X")
    return aa


# Let's dissect this function here to see what's going on and then we can make a reuseable script out of this
def translate_dna(dna: str):
    protein = ""
    for start in range(0, len(dna) - 1, 3):
        codon = dna[start : start + 3]
        aa = lookup_AA(codon)
        protein = protein + aa
    return protein


@click.command()
@click.argument("dnas")
@click.option("--ignore_stop", is_flag=True)
def main(dnas: str, ignore_stop: bool):
    for seq in dnas.split(","):
        protein = translate_dna(seq)
        print(protein.replace("_", "") if ignore_stop else protein)


if __name__ == "__main__":
    main()

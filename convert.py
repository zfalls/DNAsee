import sys
import pandas as pd

# Function to get reverse-sense DNA
def dna_complement(dna):
    # clean DNA str
    dna = dna.lower().replace(' ','').replace('\n','')
    #dna_temp = list(dna)
    #dna_temp.reverse()
    dna_temp = dna[::-1]
    comp_dna = []
    for i in dna_temp:
        if i == 't':
            comp_dna += ['a']
        elif i == 'a':
            comp_dna += ['t']
        elif i == 'g':
            comp_dna += ['c']
        elif i == 'c':
            comp_dna += ['g']
        else:
            print("Unknown nucleotide {}...exiting program".format(i))
            sys.exit()
    comp_dna = ''.join(comp_dna)
    return comp_dna

# Function to convert DNA sequence to RNA sequence
def dna2rna(dna):
    rna = []
    # clean DNA str
    dna = dna.lower().replace(' ','').replace('\n','')
    # create rna sequence with complements to dna sequence
    for i in dna:
        if i == 't':
            rna += ['u']
        elif i == 'a':
            rna += ['a']
        elif i == 'g':
            rna += ['g']
        elif i == 'c':
            rna += ['c']
        else:
            print("Unknown nucleotide {}...exiting program".format(i))
            sys.exit()
    rna = ''.join(rna)
    return rna


# Function to convert RNA sequence to DNA sequence
def rna2dna(rna):
    dna = []
    # clean DNA str
    rna = rna.lower().replace(' ','').replace('\n','')
    # create rna sequence with complements to dna sequence
    for i in rna:
        if i == 'u':
            dna += ['t']
        elif i == 'a':
            dna += ['a']
        elif i == 'g':
            dna += ['g']
        elif i == 'c':
            dna += ['c']
        else:
            print("Unknown nucleotide {}...exiting program".format(i))
            sys.exit()
    dna = ''.join(dna)
    return dna

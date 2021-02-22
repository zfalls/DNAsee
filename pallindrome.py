import pandas as pd
from convert import *
from secondary_struct import *

def check_bulge(dna,temp_i,temp_j):
    """ Check if a 1 nucelotide bulge exists at the position temp_j.

    Parameters:
        dna (str): Full RNA sequence
        temp_i (int): current left position (index of dna) of the stem-loop being analyzed
        temp_j (int): current right position (index of dna) of the stem-loop being analyzed

    Returns:
        bool: Does a 1 nucleotide bulge exist at this position
    """

    temp_j+=1
    if dna[temp_i] == 'c' and dna[temp_j] == 'g':
        return True
    elif dna[temp_i] == 'g' and dna[temp_j] == 'c':
        return True
    elif dna[temp_i] == 'a' and dna[temp_j] == 't':
        return True
    elif dna[temp_i] == 't' and dna[temp_j] == 'a':
        return True
    else:
        return False

# Function to determine the length of the pallindromic sequence
# Stem length
def pallindrome(dna,best,loop,length,pos_c,i,j,columns):
    df = pd.DataFrame()
    temp_i = i-1
    temp_j = j+1
    stem = 0
    bulge = False
    while temp_i >= 0 and temp_j < length:
        if dna[temp_i] == 'c' and dna[temp_j] == 'g':
            stem += 1
        elif dna[temp_i] == 'g' and dna[temp_j] == 'c':
            stem += 1
        elif dna[temp_i] == 'a' and dna[temp_j] == 't':
            stem += 1
        elif dna[temp_i] == 't' and dna[temp_j] == 'a':
            stem += 1
        # Check for bulge at -2 position ONLY
        elif temp_j - pos_c == 1 and check_bulge(dna,temp_i,temp_j):
            bulge = True
            temp_i -= 1
            temp_j += 2
            stem += 1
            continue
        else:
            break
        # Keep searching...
        temp_i -= 1
        temp_j += 1
    if stem < 1:
        return df
    temp_i += 1
    temp_j -= 1
    # Temp store stem-loop sequence
    temp_dna = dna[temp_i:temp_j+1]
    # Get dna complement of dna stem-loop sequence
    temp_rna = dna2rna(temp_dna)
    # Get secondary structure and energy
    ss, mfe = get_ss(temp_rna)

    # Leave commented for now...need to test more
    '''
    # Get secondary structure and energy
    # of a larger subsequence (length of 40nt)
    if bulge:
        diff = temp_j - temp_i + 1
    else:
        diff = temp_j - temp_i
    change = (40 - diff) / 2
    if temp_i - change >= 0 and temp_j + change <= len (dna):
        temp_dna_large = dna[temp_i-change:temp_j+change+1]
        ss_large, mfe = get_ss(temp_dna_large)
    elif temp_i - change < 0:
        temp_dna_large = dna[0:temp_j+1+temp_i]
        ss_large, mfe = get_ss(temp_dna_large)
    elif temp_j + change > len(dna):
        diff = len(dna) - temp_j
        temp_dna_large = dna[temp_i-diff:len(dna)+1]
        ss_large, mfe = get_ss(temp_dna_large)
    '''

    # Store all stem-loop info
    df = pd.DataFrame([[temp_dna,temp_rna,ss,best,loop,stem,bulge,pos_c,temp_i+1,temp_j+1,mfe]], columns=columns)
    return df

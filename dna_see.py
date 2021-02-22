#!/usr/bin/env python

from Bio.Seq import *
import pandas as pd
from convert import *
from pallindrome import *
from secondary_struct import *
import argparse
import RNA
import os, math

def aa_change(pos_c,prot,rna,j):
    pos = math.ceil(pos_c/3.0)
    pos = int(pos)
    # Edit the C>U
    temp_rna = list(rna)
    temp_rna[j] = 'u'
    temp_rna = ''.join(temp_rna)
    # Translate the edited sequence
    temp_prot = translate(temp_rna)
    if prot[pos-1] != temp_prot[pos-1]:
        f = "{} -> {}".format(prot[pos-1],temp_prot[pos-1])
    else:
        f = "synonymous"
    return f,pos

   
def see(seq, is_rna=False, tetraloop=False, reverse_sense=False): 
    out = seq.split("/")[-1].replace(".fasta","")
    
    # Read in input fasta
    f = open(seq,'r')
    f.readline()
    l = f.readlines()
    l = [i.strip() for i in l] 
    rna = ''

    if is_rna:
        # Input was RNA sequence
        rna = ''.join(l)
        rna = rna.lower()
        # Convert to DNA
        dna = rna2dna(rna)
    else:
        # Input was DNA
        dna = ''.join(l)
        dna = dna.lower()
        # Check sense
        if reverse_sense:
            dna = dna_complement(dna)            
        # Convert to RNA
        rna = dna2rna(dna)

    prot = translate(rna)

    # Create dataframes
    if tetraloop:
        columns=('dna','rna','secondary_struct','best_loop','loop_len','stem_len','bulge','pos_c','begin','end','mfe')
    else:
        columns=('dna','rna','secondary_struct','cc/tc','loop_len','stem_len','bulge','pos_c','begin','end','mfe')
    df_all = pd.DataFrame(columns=columns)
    df3 = pd.DataFrame(columns=columns)
    df4 = pd.DataFrame(columns=columns)
    df4_best = pd.DataFrame(columns=columns)
    df5 = pd.DataFrame(columns=columns)
    
    size = len(dna)
    i = 0
    j = i+3
    struct = 1
    if not os.path.exists(out): os.mkdir(out)
    while j < size:
        if dna[j] == 'c':
            pos_c = j+1
            if tetraloop:
                if dna[i:j+1] == 'catc' \
                        or dna[i:j+1] == 'cacc' \
                        or dna[i:j+1] == 'cctc' \
                        or dna[i:j+1] == 'cttc' \
                        or dna[i:j+1] == 'tatc':
                    loop = 4
                    best = 1
                    temp_df = pallindrome(dna,best,loop,size,pos_c,i,j,columns)
                    a,b = aa_change(pos_c,prot,rna,j)
                    temp_df['aa_change'] = a
                    temp_df['aa_pos'] = int(b)

                    if not temp_df.empty:
                        df4_best = df4_bests.append(temp_df, ignore_index=True, sort=False)
                        # Plot the 2D Structure
                        #print(temp_df.loc[0,'rna'])
                        #plot_ss(temp_df.loc[0,'rna'], temp_df.loc[0,'secondary_struct'], "{}/{}.ps".format(out,struct))
                        #struct += 1
                else:
                    # Check other tetra-loops anyway
                    best = 0
                    loop = 4
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()

                    if not temp_df.empty:
                        df4 = df4.append(temp_df, ignore_index=True, sort=False)
                    
                    # Check for tri- and penta-loops here
                    shift_i = i+1
                    loop = 3
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,shift_i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()

                    if not temp_df.empty:
                        df3 = df3.append(temp_df, ignore_index=True, sort=False)

                    shift_i = i-1
                    loop = 5
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,shift_i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()

                    if not temp_df.empty:
                        df5 = df5.append(temp_df, ignore_index=True, sort=False)
    
            else:
                if dna[j-1] == 'u' or dna[j-1] == 'c':
                    loop = 4
                    best = 1
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()
                    if not temp_df.empty:
                        df4_best = df4_best.append(temp_df, ignore_index=True, sort=False)
                        # Plot the 2D Structure
                        #print(temp_df.loc[0,'rna'])
                        #plot_ss(temp_df.loc[0,'rna'], temp_df.loc[0,'secondary_struct'], "{}/{}.ps".format(out,struct))
                        #struct += 1

                    # Check for tri- and penta-loops here
                    shift_i = i+1
                    loop = 3
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,shift_i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()
                    if not temp_df.empty:
                        df3 = df3.append(temp_df, ignore_index=True, sort=False)
                    
                    shift_i = i-1
                    loop = 5
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,shift_i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()
                    if not temp_df.empty:
                        df5 = df5.append(temp_df, ignore_index=True, sort=False)

                else:
                    # Check other tetra-loops anyway
                    best = 0
                    loop = 4
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()
                    if not temp_df.empty:
                        df4 = df4.append(temp_df, ignore_index=True, sort=False)
                        #i += 1
                        #j += 1
                        #continue

                    # Check for tri- and penta-loops here
                    shift_i = i+1
                    loop = 3
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,shift_i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()
                    if not temp_df.empty:
                        df3 = df3.append(temp_df, ignore_index=True, sort=False)
                    
                    shift_i = i-1
                    loop = 5
                    try:
                        temp_df = pallindrome(dna,best,loop,size,pos_c,shift_i,j,columns)
                        a,b = aa_change(pos_c,prot,rna,j)
                        temp_df['aa_change'] = a
                        temp_df['aa_pos'] = int(b)
                    except:
                        temp_df = pd.DataFrame()
                    if not temp_df.empty:
                        df5 = df5.append(temp_df, ignore_index=True, sort=False)
        i += 1
        j += 1
    # Sort df by stem length (descending)
    df4_best = df4_best.sort_values(by=['stem_len','mfe'], ascending=[False,True])
    df4 = df4.sort_values(by=['mfe'], ascending=True)
    df4 = df4.sort_values(by=['stem_len'], ascending=False)
    df3 = df3.sort_values(by=['stem_len','mfe'], ascending=False)
    df5 = df5.sort_values(by=['stem_len','mfe'], ascending=False)
    
    # Merge dfs 
    df_all = df_all.append([df4_best,df4,df3,df5], ignore_index=True, sort=False)
    df_all['rank'] = df_all.index+1
    df_all = df_all.astype({'aa_pos':int})
    if tetraloop:
        df_all = df_all[['rank','pos_c','begin','end','loop_len','stem_len','best_loop','bulge','dna','rna','secondary_struct','mfe','aa_change','aa_pos']]
    else:
        df_all = df_all[['rank','pos_c','begin','end','loop_len','stem_len','cc/tc','bulge','dna','rna','secondary_struct','mfe','aa_change','aa_pos']]
    
    # Final sort by MFE (testing it out)
    df_all = df_all.sort_values(by=['mfe','bulge'], ascending=True)
 
    print(df_all)
    df_all.to_csv("{}.tsv".format(out),sep='\t',float_format='%.2f')
   
    '''
    # Get MFE secondary structure and energy
    ss, mfe = get_ss(rna)
    
    # Print the 2D Structure Data
    print_ss(rna, ss, mfe, "{}.out".format(out))
    
    # Plot the 2D Structure
    plot_ss(rna, ss, "{}.ps".format(out))
    '''
    return df_all

if __name__ == '__main__':
    # Arugments
    parser = argparse.ArgumentParser(prog="RNAsee", description="Search a sequence for putative RNA editing sites by APOBEC3A/G.")
    parser.add_argument('sequence', metavar='sequence', type=str, help="FASTA file for DNA seqeunce input.")
    parser.add_argument('-r', '--reverse', action='store_true', help="Positive or negative sense for the input DNA", default=False)
    parser.add_argument('--tetraloop','-t', action="store_true", help="Search and rank based upon specific tetraloop (CUAC, CACC, CCUC, CUUC, and UAUC).", default=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s v0.2')
    parser.add_argument('--rna', action="store_true", help="The FASTA file is an RNA sequence (only A, U, G, C).")
    args=parser.parse_args()
    
    # Set variables
    seq = args.sequence
    reverse_sense = args.reverse
    is_rna = args.rna
    tetraloop = args.tetraloop
    
    df = see(seq, is_rna, tetraloop, reverse_sense)
 

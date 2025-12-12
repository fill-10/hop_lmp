##--- README ---##
# This file contains: 
# Standard HPS-Urry force field database for 20 amino acids
# It is universal, independent of your simulation.
# Loaded from FF .dat file.
# X. Luo Nov. 2025
##-------------##

import numpy as np

##-- format pre-processing ---###
# code 1 to 3 conversion
def code_1_3( seq_1):
    # dictionary tool, CONST
    SEQ_DICT  = { \
        'R':'ARG' ,'H':'HIS', 'K':'LYS', 'D':'ASP', 'E':'GLU', \
        'S':'SER' ,'T':'THR', 'N':'ASN', 'Q':'GLN', 'C':'CYS', \
        'U':'SEC' ,'G':'GLY', 'P':'PRO', 'A':'ALA', 'V':'VAL', \
        'I':'ILE' ,'L':'LEU', 'M':'MET', 'F':'PHE', 'Y':'TYR', \
        'W':'TRP' }
    # output init, empty
    seq_3 = []
    # loop to convert
    for i in seq_1:
        if i in SEQ_DICT: # ignore if not in
            seq_3.append( SEQ_DICT[i] )
    return seq_3

def code_3_1(seq_3):
    # dictionary tool, CONST
    SEQ_DICT  = { \
        'ARG':'R','HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', \
        'SER':'S','THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', \
        'SEC':'U','GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', \
        'ILE':'I','LEU':'L', 'MET':'M', 'PHE':'F', 'TYR':'Y', \
        'TRP':'W'}
    # output init, empty
    seq_1 = []
    # loop to convert
    for i in seq_3:
        if i in SEQ_DICT: # ignore if not in
            seq_1.append( SEQ_DICT[i] )
    return seq_1

##--- convert ff file to atomwise dict
def read_ff_file ( ff_file ):
    # convert txt ff file into a dict 
    # { 'atom type': np.array( [ mass, charge, sigma, lambda], dtype=float)
    with open ( ff_file , 'r' ) as f:
        ff_raw = f.readlines()
    ff_raw =  [x for x in ff_raw if not x.startswith('#')]
    atom_ff_dict = {}
    # each atom type 
    for i in ff_raw : 
        at_line = i.split() # split by space, drop \n
        at_nm = at_line[0] # extract atom name
        # mass, charge, sigma, lambda in numpy array (float)
        at_params = np.array( at_line[ 1 : ], dtype=float)
        atom_ff_dict [ at_nm ] = at_params
    return atom_ff_dict 

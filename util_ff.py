##--- README ---##
# This file contains: 
# Standard HPS-Urry force field database for 20 amino acids
# Ready to use in HooMD simulation as a utility.
# It is universal, independent of your simulation.
# Loaded from FF .dat file.
# HOWTO: "from this_file.py import *"
# X. Luo Nov. 2025
##-------------##

import numpy as np
import hoomd
from hoomd import azplugins
#import gsd, gsd.hoomd, gsd.pygsd

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

##--- build nonbonded pairwise ff from atomwise parameters
def nb_hsp_urry( atom_ff , r_cut_default = 2.0 , \
    neighbor_list = hoomd.md.nlist.Cell( buffer = 3.0 ) ):
    # 
    nb = azplugins.pair.PerturbedLennardJones( \
        default_r_cut= r_cut_default, nlist = neighbor_list )
    #pair-wise
    ff_dict = atom_ff
    for idx, keyi in enumerate ( ff_dict.keys()):
        for jdx, keyj in enumerate ( ff_dict.keys()):
            if idx<=jdx:
                lmd_i = ff_dict[keyi][ 3 ]
                lmd_j = ff_dict[keyj][ 3 ]
                sig_i = ff_dict[keyi][ 2 ]
                sig_j = ff_dict[keyj][ 2 ]
                # arithmetic
                lmd_i_j = ( lmd_i + lmd_j ) * 0.5
                sig_i_j = ( sig_i + sig_j ) * 0.5 * 0.1 # A to nm
                nb.params[ (keyi, keyj) ] = dict ( \
                    attraction_scale_factor = lmd_i_j, \
                    epsilon = 0.8368 , sigma  = sig_i_j )
                    #, r_cut = 2.0  )
    return nb

##--- build coul pairwise ff from atomwise parameters
def coul_yukawa( atom_ff,  r_cut_default = 3.0 , \
    neighbor_list = hoomd.md.nlist.Cell( buffer = 4.0 ) ):
    yukawa = hoomd.md.pair.Yukawa( \
        default_r_cut= r_cut_default, nlist = neighbor_list )
    #pair-wise
    ff_dict = atom_ff
    for idx, keyi in enumerate ( ff_dict.keys()):
        for jdx, keyj in enumerate ( ff_dict.keys()):
            if idx<=jdx:
                q_i = ff_dict[keyi][ 1 ]
                q_j = ff_dict[keyj][ 1 ]
                q_eps = q_i * q_j * 1.73136
                yukawa.params[ (keyi, keyj ) ] = dict ( \
                    epsilon = q_eps, kappa=1.0 )
                    # r_cut=3.5) 
    return yukawa

# test:
if __name__ == "__main__":
    with open("TDP-43_LC_267-414.dat", 'r') as f:
        seq_1 = ''.join( f.readlines() ).replace('\n', '')
    tdp_lc = code_1_3( seq_1) 
    print ( "total number of aa in TDP LC: ", len (tdp_lc ) )
    nl =  hoomd.md.nlist.Cell( buffer = 4.0 ) 
    ff_atomwise = read_ff_file( 'stats_module.dat' )
    ff_nb_hps = nb_hsp_urry( ff_atomwise, neighbor_list = nl )
    print  (ff_nb_hps)
    print  ( len (ff_nb_hps.params) )
    ff_q_yukawa = coul_yukawa( ff_atomwise, neighbor_list = nl )
    print ( ff_q_yukawa )
    print ( ff_q_yukawa.params[ ('HIS', 'HIS') ]  ) 

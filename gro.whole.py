from class_data import data
if __name__ == '__main__':
    import time as timer
    import numpy as np
    import pandas as pd 
    from pbc_dist import pbc_vec
    ##--- output general settings ---
    fn_prefix = '../'
    ##--- timer start ---
    start = timer.perf_counter()
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0 #do not clear L_atoms

    ##--- read in gro traj ---
    grofilename = '../trj_P4_nojump.gro'
    d1.read_all_gro(grofilename )
    #d2 = data()
    #d2.save_mem=0 #do not clear L_atoms
    #grofilename = '/global/scratch/users/xluo2/test_charmm/Ndc10-Nte10/TFA_urea_ld/9layers/trj_P3.gro'
    #d2.read_all_gro(grofilename )
    #d1.allframes += d2.allframes[1:]
    #del d2
    #print('last timestep read: ', d1.allframes[-1].time)
    print( 'number of frames: ', len(d1.allframes) )
    ##--- give newer molid ---
    Nchain = 12 * 4
    Napc = 685 # atoms per chain
    # select Y direction: a very rough selection.
    # be carefull! The selection is based on visuallization.
    selrowy  =  (  d1.allframes[0].L_atom['yu']  > d1.allframes[0].deltaY *0.75 ) \
      &  (   ( d1.allframes[0].L_atom['res'] == 'NDC' ) \
           | ( d1.allframes[0].L_atom['res'] == 'NTE' ) \
         )
    # select Z direction
    selrowz  =  (  d1.allframes[0].L_atom['zu']  < d1.allframes[0].deltaZ *0.2 ) \
      &  (   ( d1.allframes[0].L_atom['res'] == 'NDC' ) \
           | ( d1.allframes[0].L_atom['res'] == 'NTE' ) \
         )

    for frame in d1.allframes:
        ##--- delete UREA and TFA to save time ---
        #    frame.L_atom = frame.L_atom[ \
        #    ( frame.L_atom['res'] == 'NDC' )  | \
        #    ( frame.L_atom['res'] == 'NTE' )  ]
        #
        ##-- wrap discontinuous molecules
        # move Y
        frame.L_atom.loc[ selrowy, 'yu'] -= frame.deltaY
        # move Z
        frame.L_atom.loc[ selrowz, 'zu'] += frame.deltaZ

        ##--- unit conversion ---
        """
        frame.L_atom['xu'] *= 10
        frame.L_atom['yu'] *= 10
        frame.L_atom['zu'] *= 10
        frame.xlo *= 10
        frame.xhi *= 10
        frame.ylo *= 10
        frame.yhi *= 10
        frame.zlo *= 10
        frame.zhi *= 10
        frame.deltaX *= 10
        frame.deltaY *= 10
        frame.deltaZ *= 10
        """
        ##----------------------- 
    # export
    d1.export_all_gro(fn_prefix+'trj_P4.gro')
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


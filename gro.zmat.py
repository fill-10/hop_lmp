from class_data import data
from class_oneframe import oneframe

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
    grofilename = '../test_P2.gro'
    #grofilename = './P2_NdcNoH_original.gro'
    d1.read_all_gro(grofilename )
    #
    Natom = 140
    Nframe =10002
    print('last timestep read: ', d1.allframes[-1].time)
    print( 'number of frames: ', len(d1.allframes) )
    
    d1.gen_zmat( is_sqrt= True )

    for frame in d1.allframes:
        print(len(frame.L_zmat) )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


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
    grofilename = '../trj_P3_pep_urea.gro'
    d1.read_all_gro(grofilename )
    d2 = data()
    d2.save_mem=0 #do not clear L_atoms
    grofilename = '../trj_P4_pep_urea.gro'
    d2.read_all_gro(grofilename )
    d1.allframes += d2.allframes[1:]
    del d2 # save memory
    print('last timestep read: ', d1.allframes[-1].time)
    print( 'number of frames: ', len(d1.allframes) )
    for i in range(0, len(d1.allframes)):
        f = open('../SASA/trj_P3P4_pep_urea'+str(i)+'.gro', 'w')
        d1.allframes[i].export_gro( f, d1.allframes[i].L_atom )
        f.close()
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


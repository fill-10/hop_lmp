from class_data import data
import time as timer
import numpy as np
import pandas as pd 

def min_ion_dist(d1, ion1_type, ion2_type, bin_step, output_fn):
    ##--- calcuate histograms ---
    dist = [] # all distances
    # find distance to the nearest NK atom for each CLA
    for i in range( 0, len(d1.allframes) ):
        # generate two atom list
        cr1 = d1.allframes[i].L_atom['type'] == ion1_type
        sel1 = d1.allframes[i].L_atom[ cr1 ].reset_index(drop = True )
        cr2 = d1.allframes[i].L_atom['type'] == ion2_type
        sel2 = d1.allframes[i].L_atom[ cr2 ].reset_index(drop = True)
        # for each CLA, match with all  NK
        CL_NK_min = []
        for j in range(0, sel1.shape[0] ):
            # make a data frame for one CLA
            sel3 = pd.DataFrame( \
            np.tile( sel1.iloc[j,:] , (sel2.shape[0], 1) ), \
                columns = sel1.columns )
            CL_NK_min.append( \
                min( d1.allframes[i].bond_w(sel3, sel2) )  )
                # min of a list
        dist += CL_NK_min
    bond_bins = np.arange( int(min(dist)*10)/10, int(max(dist)*10)/10+bin_step*2, bin_step)
    col_hist, col_bins = np.histogram( dist , bins=bond_bins, density=True )
    col_bins += bin_step/2
    np.savetxt( output_fn, np.transpose([col_bins[:-1], col_hist]), fmt=['%f', '%f'], header='r/nm, P(r/nm)' )


if __name__ == '__main__':
    ##--- output general settings ---
    fn_prefix = '../'
    ##--- timer start ---
    start = timer.perf_counter()
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0 #do not clear L_atoms
    ##--- read in gro traj ---
    grofilename = '../nptPR3_whole.gro'
    d1.read_all_gro(grofilename )
    d2 = data()
    d2.save_mem=0
    d2.read_all_gro('../nptPR2_whole.gro')
    d1.allframes += d2.allframes
    del d2
    print(len( d1.allframes) )
    print('last timestep read: ', d1.allframes[-1].time)
    d1.wrapall_L()
    min_ion_dist( d1, 'CLA', 'NK', 0.1, '../dist_CLA_NK_1.dat' )
    min_ion_dist( d1, 'CLA', 'N', 0.1, '../dist_CLA_Nbb_1.dat' )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)



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
    grofilename = '../trj_P4.gro'
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
    Nchain = 12 * 6
    Napc = 685 # atoms per chain
    for frame in d1.allframes:
        for i in range(0,Nchain):
            frame.L_atom.loc[Napc * i : (Napc)*(i+1)-1 ,'mol'] \
            = i+1
            # .loc[:684] includes index=684, [:684] does not.
        ##--- delete UREA and TFA to save time ---
            frame.L_atom = frame.L_atom[ \
            ( frame.L_atom['res'] == 'NDC' )  | \
            ( frame.L_atom['res'] == 'NTE' )  ]
        ##--- unit conversion ---
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
        ##----------------------- 

    d1.wrapall_L()
    #d1.export_all_lmptrj(fn_prefix+'trj_P4.lammpstrj')
    sel1_kw = [ [1, 24, 25, 48, 49, 72], \
                [1, 2, 681, 682, 683],\
                "sel[( sel.index == 370 )]"  ]
    sel2_kw = [ [1, 24, 25, 48, 49, 72], \
                [1, 2, 681, 682, 683],\
                "sel[( sel.index == 649 )]"  ]
    col_hist_dist, col_rb = d1.bond_stat( \
    sel1_kw, sel2_kw, bond_bins=np.arange(6.0, 30.0, 0.1) )
    ##--- save results ---
    np.savetxt(  fn_prefix+'dist_NTEe2e_1st.dat', \
                 np.transpose( [col_rb, col_hist_dist] ),\
                 fmt=['%f', '%f'], header='rb  P(rb)' )
    ##--- print E(r) ---
    print( 'Expectation of r for 1st: ', \
            np.sum( col_rb*col_hist_dist ) )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


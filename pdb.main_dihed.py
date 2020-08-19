from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    pdbfilename = '../testCtSt.pdb' # hopping calculation
    ##--- output general settings ---
    import numpy as np
    fn_prefix = '../AmC2_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0

    ##--- read in lmp traj output ---
    #AN_gen_kw =  [ range(11,11+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    #CT_gen_kw =  [ range(1,1+10), 2, 40, 1, [30,1171], 15, "sel[ ((sel.index%30<=9 ) & (sel.index%30>=5)) | ( (sel.index%30>=15 ) & (sel.index%30<=24)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename )#, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
    #d1.wrapall_L()

    ##--- calcuate histograms ---
    sel1_kw = [ range(11,11+400), [], "sel[ sel.index%15==7 ]" ]
    sel2_kw = [ range(11,11+400), [], "sel[ sel.index%15==2 ]" ]
    sel3_kw = [ range(11,11+400), [], "sel[ sel.index%15==0 ]" ]
    sel4_kw = [ range(11,11+400), [], "sel[ sel.index%15==8 ]" ]

    #col_hist_bond, col_dist = d1.bond_stat(sel1_kw, sel2_kw, np.arange(1,2, 0.1))
    #print(col_hist_bond, '\n*\n', col_dist)

    #col_hist_angle, col_dist = d1.angle_stat(sel1_kw, sel2_kw, sel3_kw)
    #print(len(col_hist_angle), '\n*\n', len(col_dist) )

    col_hist_dihed, col_dist = d1.dihed_stat(sel1_kw, sel2_kw, sel3_kw, sel4_kw)

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

    np.savetxt( fn_prefix+'dihed.dat', np.transpose([col_dist, col_hist_dihed]), fmt=['%f', '%f'], header='dihedral, P(phi)' )


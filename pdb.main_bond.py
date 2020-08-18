from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    lmpfilename = '/home/xubo/LiuSum/tf2n_luo/C2VIm/500/ions_every100ps_10-50ns.lammpstrj'
    ##--- output general settings ---
    import numpy as np
    fn_prefix = '../Tf2C2_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0

    ##--- read in pdb traj output ---
    #AN_gen_kw =  [ range(11,11+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    #CT_gen_kw =  [ range(1,1+10), 2, 40, 1, [30,1171], 15, "sel[ ((sel.index%30<=9 ) & (sel.index%30>=5)) | ( (sel.index%30>=15 ) & (sel.index%30<=24)) ]", 'type'   ]
    d1.read_all_lmp(lmpfilename)
    print('last timestep read: ', d1.allframes[-1].time)
    d1.unwrapall_L()

    ##--- calcuate histograms ---
    sel1_kw = [ range(401,401+10), [20, 781], "sel[ sel.index%20==8 ]" ]
    sel2_kw = [ range(401,401+10), [20, 781], "sel[ sel.index%20==0 ]" ]
    sel3_kw = [ range(811,811+10), [], "sel[:]" ]

    col_hist_bond, col_dist = d1.bond_stat(sel2_kw, sel3_kw, np.arange(0,10, 0.1))
    print(col_hist_bond, '\n', col_dist)

    col_hist_angle, col_ang = d1.angle_stat(sel1_kw, sel2_kw, sel3_kw)
    print(col_hist_angle, '\n', col_dist )

    #col_hist_dihed, col_dist = d1.dihed_stat(sel1_kw, sel2_kw, sel3_kw, sel4_kw)

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

    np.savetxt( fn_prefix+'bond.dat', np.transpose([col_dist, col_hist_bond]), fmt=['%f', '%f'], header='bond, P(r/A)' )
    np.savetxt( fn_prefix+'angle.dat', np.transpose([col_ang, col_hist_angle]), fmt=['%f', '%f'], header='angle, P(Theta/degree)' )


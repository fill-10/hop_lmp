from class_data import data
if __name__ == '__main__':
    import time as timer
    pdbfilename = '../nongauss_every10ps_00-50ns.pdb' # vanhove_s calculation
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = '../AmC2_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(11,11+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(1,1+10), 2, 40, 1, [30,1171], 15, "sel[ ((sel.index%30<=9 ) & (sel.index%30>=5)) | ( (sel.index%30>=15 ) & (sel.index%30<=24)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
    
    ##--- van hove self ---
    #d1.unwrapall_AN()  # no need to unwrap if read from fix file and uxyz not deleted
    #print('unwrapped data')
    dist_col, vanhove_s_col, fpi_r2_vanhove_s_col = d1.fpi_r2_vanhove_s_AN_avg(35, 25.0, 0.1, 5) # t*=350 ps, so interval*=35, skip = 5 frames (default skip = 0)
    print('van hove self calculated')
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)
    
    ##--- save vanhove func
    np.savetxt(    fn_prefix+'vanhove_s.dat', np.transpose( np.vstack( (dist_col, vanhove_s_col) )), fmt=['%f', '%f'], header='dist, vanhove_s'      )
    np.savetxt(    fn_prefix+'4pi_r2_vanhove_s.dat', np.transpose( np.vstack( (dist_col, fpi_r2_vanhove_s_col) )), fmt=['%f', '%f'], header='dist, 4pi*r^2vanhove_s'      )

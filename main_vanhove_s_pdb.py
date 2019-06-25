from class_data import data
if __name__ == '__main__':
    import time as timer

    pdbfilename = '../vanhove580_2200.pdb' # vanhove_s calculation
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = 'BrC5_pdb_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
    
    ##--- van hove self ---
    d1.unwrapall_AN()  # no need to unwrap if read from fix file and uxyz not deleted
    print('unwrapped data')
    dist_col, vanhove_s_col, fpi_r2_vanhove_s_col = d1.fpi_r2_vanhove_s_AN_avg(1, 25.0, 0.1) # t*=2200 interval in pdb is 550 so interval*=4
    print('van hove self calculated')
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)
    
    ##--- save vanhove func
    np.savetxt(    fn_prefix+'vanhove_s.dat', np.transpose( np.vstack( (dist_col, vanhove_s_col) )), fmt=['%f', '%f'], header='dist, vanhove_s'      )
    np.savetxt(    fn_prefix+'4pi_r2_vanhove_s.dat', np.transpose( np.vstack( (dist_col, fpi_r2_vanhove_s_col) )), fmt=['%f', '%f'], header='dist, 4pi*r^2vanhove_s'      )

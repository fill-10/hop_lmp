from class_data import data
if __name__ == '__main__':
    import time as timer
    pdbfilename = '../nongauss_every100ps.pdb' # vanhove_s calculation
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = '../Tf2C2_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
    
    ##--- wrap ---
    d1.wrapall_L()  # no need to unwrap if read from fix file and uxyz not deleted
    #print('unwrapped data')
    ##--- van hove self ---
    dist_col, vanhove_d_col, vanhove_d_o_4pir2_col = d1.vanhove_d_AN_o_4pir2_avg(13, 25.0, 0.1) # t*=1200 ps, so interval*=12
    print('van hove dist calculated')
    ##--- timer stop ---
    stop = timer.perf_counter()
    
    ##--- save vanhove func
    np.savetxt(    fn_prefix+'vanhove_d.dat', np.transpose( np.vstack( (dist_col, vanhove_d_col) )), fmt=['%f', '%f'], header='dist, vanhove_d'      )
    np.savetxt(    fn_prefix+'vanhove_d_o_4pir2.dat', np.transpose( np.vstack( (dist_col, vanhove_d_o_4pir2_col ) )), fmt=['%f', '%f'], header='dist, vanhove_d/4pi*r^2'      )
    print('time used in sec: ', stop-start)

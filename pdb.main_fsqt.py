from class_data import data
if __name__ == '__main__':
    import time as timer
    pdbfilename = '../nongauss_every100ps.pdb'
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
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw )#, CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN) 
    print('last timestep read: ', d1.allframes[-1].time)
 
    ##--- self part of intermediate scattering function ---
    time_col, fsqt_col = d1.fsqt_AN_avg(0.76, 100, 3000) # q*=0.76A
    print('fsqt calculated')
    
    ##--- timer stop ---
    stop = timer.perf_counter()

    ##--- save non gaussian parameter
    np.savetxt(    fn_prefix+'fsqt_forward_3000.dat', np.transpose( [time_col, fsqt_col] ), fmt=['%d', '%f'], header='time/ps, fsqt'      )
    print('time used in sec: ', stop-start)

from class_data import data
if __name__ == '__main__':
    import time as timer
    pdbfilename = '../nvtlong_every100ps_0-300ns.pdb' # rdf calculation
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = 'BrC5_0-50ns_lmptrj_'
    
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
 
    ##--- non gaussian parameter ---
    ##--- unwrap anion coordinates --- # may not necessary for some cases
    d1.unwrapall_AN()  # no need to unwrap if read from fix file and uxyz not deleted
    print('unwrapped data')
    timestep_col, nongauss_col = d1.nongauss_AN_avg(1, 100,1000) # start_interval = 1 frame,  100ps/frame, maxattemp =500 by default
    print('non gauss calculated')
    
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

    ##--- save non gaussian parameter
    np.savetxt(    fn_prefix+'nongauss.dat', np.transpose( [timestep_col, nongauss_col] ), fmt=['%d', '%f'], header='time/ps, nongauss'      )

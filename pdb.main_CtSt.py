from class_data import data

if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    pdbfilename = '../nvtprod_every10ps_10-20ns.pdb'
    pdbfilename = '../testCtSt.pdb'
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = '../AmC2_'
 
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1 # save atoms and ions

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(11,11+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(1,1+10), 2, 40, 1, [30,1171], 15, "sel[ ((sel.index%30<=9 ) & (sel.index%30>=5)) | ( (sel.index%30>=15 ) & (sel.index%30<=24)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    print('last timestep read: ', d1.allframes[-1].time)
    d1.wrapall_L()
    
    t1 = timer.perf_counter()
    ##--- ht ---
    d1.ht_gen(9.0)
    t2 = timer.perf_counter()
    ##--- Ct ---
    time_col, Ct_col = d1.Ct( 10, 5 )
    t3 = timer.perf_counter()
    np.savetxt(    fn_prefix+'Ct.dat', np.transpose( [time_col, Ct_col] ), fmt=['%d', '%f'], header='time/ps, Ct'      )
    ##--- St ---
    time_col, St_col = d1.St( 10, 5 )
    np.savetxt(    fn_prefix+'St.dat', np.transpose( [time_col, St_col] ), fmt=['%d', '%f'], header='time/ps, St'      )
    t4 = timer.perf_counter()
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

    print('\nCOM and wrap:',t1-start, '\nht: ', t2-t1, '\nCt: ', t3-t2, '\nSt: ', t4-t3)

from class_data import data

if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    pdbfilename = '../nvtprod_every10ps_10-20ns.pdb'

    ##--- output general settings ---=
    import numpy as np
    fn_prefix = './'
 
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
    d1.ht_gen(9.0)
    ##--- ht ---
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

    np.savetxt(    fn_prefix+'ht1.dat', d1.L_all_ht, fmt='%d')

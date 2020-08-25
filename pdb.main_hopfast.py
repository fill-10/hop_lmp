from class_data import data
if __name__ == '__main__':
    import time as timer
    pdbfilename = '../nongauss_every100ps.pdb'
    ##--- output general settings ---
    import numpy as np
    fn_prefix = '../Tf2C2_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(401,401+10), 2, 40, 1, [20, 781], 8, "sel[ ((sel.index%20<=5 ) & (sel.index%20>=1)) | ( (sel.index%20>=9 ) & (sel.index%20<=11)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    print('last timestep read: ', d1.allframes[-1].L_atom)
    d1.wrapall_L()

    ##--- calcuate histograms ---
    norm_hist_atom, norm_hist_mol = d1.find_asso_AN_CT(7.8) # cut=7.8, skip=0
    fast_percentage = d1.find_AN_fast(11, 5.9, skip=0)  # interval_star= 11= 1100ps, r*=5.9A
    norm_hist_hop_type = d1.hopfast_AN(11) # t*=1100ps=11intervals, skip=0

    ##--- save hopping types
    #np.savetxt(    fn_prefix+'hist_asso_atom.dat', np.transpose( [norm_hist_atom[1][:-1], norm_hist_atom[0][:]] ), fmt=['%d', '%f'], header='n  P(n)'      )
    #np.savetxt(    fn_prefix+'hist_asso_mol.dat' , np.transpose( [norm_hist_mol[1][:-1] , norm_hist_mol[0][:] ] ), fmt=['%d', '%f'], header='N  P(N)'      )
    np.savetxt(    fn_prefix+'hop_fast.dat'  , np.transpose( [norm_hist_hop_type[1][:-1], norm_hist_hop_type[0][:] ]  ), fmt=['%d', '%f'] , header='N  P(N)'  )

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

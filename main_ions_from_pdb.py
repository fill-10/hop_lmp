from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    pdbfilename = '../nvtlong_every10ns_0-300ns.pdb'
    #pdbfilename = '../nvtlong_every100ps_0-50ns.pdb' # rdf calculation
    pdbfilename = '../nvtlong_every10ps_0-50ns.pdb' # hopping calculation
    ##--- output general settings ---
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
    
    ##--- save lammpstrj of ions to calculate rdf---
    """
    d1.export_ions_lmptrj('ions_every100ps_0-50ns.lammpstrj')
    print('ions written in lammpstrj')
    """

    ##--- calcuate histograms ---
    #"""
    norm_hist_atom, norm_hist_mol = d1.find_asso_AN_CT(6.3)
    norm_hist_hop_type = d1.hoppingtype_AN()

    ##--- save hopping types
    np.savetxt(    fn_prefix+'hist_asso_atom.dat', np.transpose( [norm_hist_atom[1][:-1], norm_hist_atom[0][:]] ), fmt=['%d', '%f'], header='n  P(n)'      )
    np.savetxt(    fn_prefix+'hist_asso_mol.dat' , np.transpose( [norm_hist_mol[1][:-1] , norm_hist_mol[0][:] ] ), fmt=['%d', '%f'], header='N  P(N)'      )
    np.savetxt(    fn_prefix+'hist_hopping.dat'  , np.transpose( [norm_hist_hop_type[1][:-1], norm_hist_hop_type[0][:] ]  ), fmt=['%d', '%f'] , header='N  P(N)'  )
    #"""
    ##--- non gaussian parameter ---
    ##--- unwrap anion coordinates --- # may not necessary for some cases
    #d1.unwrapall_AN()
    #print('unwrapped AN coordinates')
    #timestep_col, nongauss_col = d1.nongauss_AN_avg(1, 100) # 100ps/frame
    #print('non gauss calculated')

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)
    
    ##--- save non gaussian parameter
    """
    np.savetxt(    fn_prefix+'nongauss.dat', np.transpose( [timestep_col, nongauss_col] ), fmt=['%d', '%f'], header='timestep, nongauss'      )
    """

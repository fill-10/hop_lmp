from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    filename = ''
    anionfixfile = '../400K_corrected_lmp/com_AN_every10ps_0-10ns.dat'
    cationfixfile = '../400K_corrected_lmp/com_CT_every10ps_0-10ns.dat'
    anionfixfile = '../com_AN_every10ps_0-10ns.dat'
    cationfixfile = '../com_CT_every10ps_0-10ns.dat'

    ##--- timer start ---
    start = timer.perf_counter()
    
    d1 = data(filename)
    d1.AN_CT_readfix(  anionfixfile, 1, cationfixfile, 40, [ [-0.02,58.5387],[-0.02,58.5387],[0.,58.5187]  ]  )
    d1.export_ions_lmptrj('ions.lammpstrj')
    
    ##--- calcuate histograms ---
    norm_hist_atom, norm_hist_mol = d1.find_asso_AN_CT(7.8)
    norm_hist_hop_type = d1.hoppingtype_AN()
    
    ##--- timer stop ---
    stop = timer.perf_counter()
    print(stop-start)
    
    ##--- output ---
    import numpy as np
    fn_prefix = 'C4_0-20ns'
    np.savetxt(    fn_prefix+'hist_asso_atom.dat', np.transpose( [norm_hist_atom[1][:-1], norm_hist_atom[0][:]] ), fmt=['%d', '%f'], header='n  P(n)'      )
    np.savetxt(    fn_prefix+'hist_asso_mol.dat' , np.transpose( [norm_hist_mol[1][:-1] , norm_hist_mol[0][:] ] ), fmt=['%d', '%f'], header='N  P(N)'      )
    np.savetxt(    fn_prefix+'hist_hopping.dat'  , np.transpose( [norm_hist_hop_type[1][:-1], norm_hist_hop_type[0][:] ]  ), fmt=['%d', '%f'] , header='type  P(type)'  )

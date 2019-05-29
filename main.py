from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    lmpfilename = '../400K_corrected_lmp_B/VImC4_every100ps_0-50ns.lammpstrj'
    anionfixfile = '../400K_corrected_lmp_B/com_AN_every100ps_0-50ns.dat'
    cationfixfile = '../400K_corrected_lmp_B/com_CT_every100ps_0-50ns.dat'

    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    
    ##--- read in ions from lammps 'fix' output ---
    d1.AN_CT_readfix(  anionfixfile, 1, cationfixfile, 40, [ [0.0,55.8983],[0.02,55.9183],[0.01,55.9083]  ], False  )
    d1.export_ions_lmptrj('ions_lmpfix.lammpstrj')
    print('ions written in lammpstrj') 

    ##--- calcuate histograms ---
    #norm_hist_atom, norm_hist_mol = d1.find_asso_AN_CT(7.8)
    #norm_hist_hop_type = d1.hoppingtype_AN()

    ##--- non gaussian parameter ---
    ##--- unwrap anion coordinates --- # may not necessary for some cases
    d1.unwrapall_AN()
    print('unwrapped data')
    timestep_col, nongauss_col = d1.nongauss_AN_avg(1, 100) # 100ps/frame
    print('non gauss calculated')

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)
    
    ##--- output ---
    import numpy as np
    fn_prefix = 'C2_0-50ns_lmpfix_'

    ##--- save hopping types
    """
    np.savetxt(    fn_prefix+'hist_asso_atom.dat', np.transpose( [norm_hist_atom[1][:-1], norm_hist_atom[0][:]] ), fmt=['%d', '%f'], header='n  P(n)'      )
    np.savetxt(    fn_prefix+'hist_asso_mol.dat' , np.transpose( [norm_hist_mol[1][:-1] , norm_hist_mol[0][:] ] ), fmt=['%d', '%f'], header='N  P(N)'      )
    np.savetxt(    fn_prefix+'hist_hopping.dat'  , np.transpose( [norm_hist_hop_type[1][:-1], norm_hist_hop_type[0][:] ]  ), fmt=['%d', '%f'] , header='N  P(N)'  )
    """

    ##--- save non gaussian parameter
    np.savetxt(    fn_prefix+'nongauss.dat', np.transpose( [timestep_col, nongauss_col] ), fmt=['%d', '%f'], header='timestep, nongauss'      )

from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    lmpfilename = '../400K_corrected_lmp_B/VImC4_every100ps_0-50ns.lammpstrj'
    anionfixfile = '../400K_corrected_lmp_B/com_AN_every100ps_0-50ns.dat'
    cationfixfile = '../400K_corrected_lmp_B/com_CT_every100ps_0-50ns.dat'
    #lmpfilename = './test10.lammpstrj'
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d2 = data()
    
    ##--- read in ions from lammps 'fix' output ---
    #d1.AN_CT_readfix(  anionfixfile, 1, cationfixfile, 40, [ [-0.02,58.5387],[-0.02,58.5387],[0.,58.5187]  ], False  )
    #d1.export_ions_lmptrj('ions.lammpstrj')
    #print('ions written in lammpstrj') 

    ##--- read in lammps traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[:]" ]
    CT_gen_kw =  [ range(411,411+400), 2, 1, 1, [], 8 , "sel[:]"]
    #CT_gen_kw =  [ range(401,401+10), 2, 40, 1, [20, 781], 8,  "sel[ ( (sel.index%20<=5)&(sel.index%20>=1) ) | ( (sel.index%20>=9)&(sel.index%20<=11)  )  ]" ]
    d2.read_all_lmp(lmpfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    # correct mol 
    NDeg_Poly_CT = 40
    NIon_permon = 1
    Nmol_CT = 10
    for counter, frame in enumerate(d2.allframes):
        for j0 in range(0, Nmol_CT):
            for j1 in range(j0*NDeg_Poly_CT*NIon_permon, (j0+1)*NDeg_Poly_CT*NIon_permon):
                frame.L_CT.loc[j1, 'mol'] = j0+1+frame.L_AN.shape[0]
    d2.export_ions_lmptrj('ions_lmptrj.lammpstrj')
    print('ions written in lammpstrj') 
    
    ##--- calcuate histograms ---
    #norm_hist_atom, norm_hist_mol = d2.find_asso_AN_CT(7.8)
    #norm_hist_hop_type = d2.hoppingtype_AN()

    ##--- non gaussian parameter ---
    ##--- unwrap anion coordinates --- # may not necessary for some cases
    d2.unwrapall_AN()
    print('unwrapped AN coordinates')
    timestep_col, nongauss_col = d2.nongauss_AN_avg(1, 100) # 100ps/frame
    print('non gauss calculated')

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)
    
    ##--- output ---
    import numpy as np
    fn_prefix = 'C2_0-50ns_lmptrj_'

    ##--- save hopping types
    """
    np.savetxt(    fn_prefix+'hist_asso_atom.dat', np.transpose( [norm_hist_atom[1][:-1], norm_hist_atom[0][:]] ), fmt=['%d', '%f'], header='n  P(n)'      )
    np.savetxt(    fn_prefix+'hist_asso_mol.dat' , np.transpose( [norm_hist_mol[1][:-1] , norm_hist_mol[0][:] ] ), fmt=['%d', '%f'], header='N  P(N)'      )
    np.savetxt(    fn_prefix+'hist_hopping.dat'  , np.transpose( [norm_hist_hop_type[1][:-1], norm_hist_hop_type[0][:] ]  ), fmt=['%d', '%f'] , header='N  P(N)'  )
    """

    ##--- save non gaussian parameter
    np.savetxt(    fn_prefix+'nongauss.dat', np.transpose( [timestep_col, nongauss_col] ), fmt=['%d', '%f'], header='timestep, nongauss'      )

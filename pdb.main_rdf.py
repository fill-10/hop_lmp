from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    #pdbfilename = '../nvtlong_every100ps_0-50ns.pdb' # rdf calculation
    pdbfilename = '../nvtprod_every1ns_00-50ns.pdb' # hopping calculation
    ##--- output general settings ---
    import numpy as np
    fn_prefix = 'ImC5_00-50ns_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
    for frame in d1.allframes: #adjust id and mol
        frame.L_AN['id'] += max(frame.L_atom['id'])
        frame.L_AN['mol'] += max(frame.L_atom['mol'])
        frame.L_CT['id'] += max(frame.L_atom['id'])
        frame.L_CT['mol'] += max(frame.L_atom['mol'])
    d1.wrapall_L()
    ##--- save lammpstrj of ions to calculate rdf---
    d1.export_all_lmptrj('ions_every1ns_0-50ns.lammpstrj')
    d1.export_all_pdb('ions_every1ns_0-50ns.pdb')
    print('ions written in lammpstrj')

    ##--- calcuate histograms ---
    """
    norm_hist_atom, norm_hist_mol = d1.find_asso_AN_CT(6.3)
    norm_hist_hop_type = d1.hoppingtype_AN()

    ##--- save hopping types
    np.savetxt(    fn_prefix+'hist_asso_atom.dat', np.transpose( [norm_hist_atom[1][:-1], norm_hist_atom[0][:]] ), fmt=['%d', '%f'], header='n  P(n)'      )
    np.savetxt(    fn_prefix+'hist_asso_mol.dat' , np.transpose( [norm_hist_mol[1][:-1] , norm_hist_mol[0][:] ] ), fmt=['%d', '%f'], header='N  P(N)'      )
    np.savetxt(    fn_prefix+'hist_hopping.dat'  , np.transpose( [norm_hist_hop_type[1][:-1], norm_hist_hop_type[0][:] ]  ), fmt=['%d', '%f'] , header='N  P(N)'  )

    """

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

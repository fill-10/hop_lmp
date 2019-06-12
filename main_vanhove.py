from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    lmpfilename = '../400K_corrected_lmp/VImC4_every100ps_0-50ns.lammpstrj'
    anionfixfile = '../400K_corrected_lmp/com_AN_every100ps_0-50ns.dat'
    cationfixfile = '../400K_corrected_lmp/com_CT_every100ps_0-50ns.dat'

    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d2 = data()
    
    ##--- read in ions from lammps 'fix' output ---
    d1.AN_CT_readfix(  anionfixfile, 1, cationfixfile, 40, [ [0.0,55.8983],[0.02,55.9183],[0.01,55.9083]  ], False  )  # an_file, AN_ionspermol, ct_file, CT_ionspermol, box, dropuxyz = False
    #d1.export_ions_lmptrj('ions_lmpfix.lammpstrj')
    #print('ions written in lammpstrj') 
    print('file:\n' + anionfixfile +'\n' + cationfixfile +'\n' +'were loaded.') 
    
    ##--- append more time frames ---
    an_fixfile = '../400K_corrected_lmp/continue/com_AN_every100ps_50-80ns.dat'
    ct_fixfile = '../400K_corrected_lmp/continue/com_CT_every100ps_50-80ns.dat'
    d2 = data()
    d2.AN_CT_readfix(  an_fixfile, 1, ct_fixfile, 40, [ [0.0,55.8983],[0.02,55.9183],[0.01,55.9083]  ], False  )
    d1.allframes += d2.allframes[1:]
    del d2
    print('file:\n' + an_fixfile +'\n' + ct_fixfile +'\n' +'were loaded.') 
    d2 = data()
    an_fixfile = '../400K_corrected_lmp/continue2/com_AN_every100ps_80-120ns.dat'
    ct_fixfile = '../400K_corrected_lmp/continue2/com_CT_every100ps_80-120ns.dat'
    d2.AN_CT_readfix(  an_fixfile, 1, ct_fixfile, 40, [ [0.0,55.8983],[0.02,55.9183],[0.01,55.9083]  ], False  )
    d1.allframes += d2.allframes[1:]
    del d2
    print('file:\n' + an_fixfile +'\n' + ct_fixfile +'\n' +'were loaded.') 

    ##--- van hove self ---
    dist_col, fpi_r2_vanhove_s_col = d1.fpi_r2_vanhove_s_AN_avg(1000, 25.0, 0.1)
    print('van hove self calculated')
    ##--- van hove distinct ---
    dist_col2, vanhove_d_col = d1.vanhove_d_AN_avg(1000, 25.0, 0.1)
    print('van hove distinct calculated')
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)
    
    ##--- output ---
    import numpy as np
    fn_prefix = 'C2_lmpfix_'

    ##--- save vanhove func
    np.savetxt(    fn_prefix+'4pi_r2_vanhove_s.dat', np.transpose( np.vstack( (dist_col, fpi_r2_vanhove_s_col) )), fmt=['%f', '%f'], header='dist, 4pi*r^2vanhove_s'      )
    np.savetxt(    fn_prefix+'vanhove_d.dat', np.transpose( np.vstack( (dist_col2, vanhove_d_col) )), fmt=['%f', '%f'], header='dist, vanhove_d'      )

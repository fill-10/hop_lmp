from class_data import data
if __name__ == '__main__':
    import time as timer
    import numpy as np
    import pandas as pd 
    from pbc_dist import pbc_vec
    ##--- read-in the whole file ---
    grofilename = '/Users/xuboluo/Documents/test_charmm/Ndc10-Nte10/H-Ndc-Nte-NH2/rest/T300/resttrj.gro'
    ##--- output general settings ---
    fn_prefix = '../rest_T300_'
    ##--- timer start ---
    start = timer.perf_counter()
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0 #do not clear L_atoms

    ##--- read in gro traj ---
    d1.read_all_gro(grofilename )
    print('last timestep read: ', d1.allframes[-1].time)
    #for frame in d1.allframes:
    #    print(frame.L_atom)
    sel1_kw = [[], [], "  (self.L_atom['type'] == 'CA5')  \
                          & (self.L_atom['mol'] >0  )"  ]
    sel2_kw = [[], [], "  (self.L_atom['type'] == 'N') \
                          & (self.L_atom['mol'] >0  )"]
    sel3_kw = [[], [], "  (self.L_atom['type'] == 'CA')  \
                          & (self.L_atom['mol'] >0  )"]
    sel4_kw = [[], [], "  (self.L_atom['type'] == 'C')  \
                          & (self.L_atom['mol'] >0   )"]
    
    """
    print(  d1.allframes[-1].L_atom[ \
                (d1.allframes[-1].L_atom['type'] =='N' ) \
                & (d1.allframes[-1].L_atom['mol'] <10 ) \
                                   ]  \
         )
    sel1_kw = [[], [], " self.L_atom['mol'] <10"]
    sel2_kw = [[], [], " self.L_atom['mol'] <10"]
    sel3_kw = [[], [], " self.L_atom['mol'] <10"]
    sel4_kw = [[], [], " self.L_atom['mol'] >1"]
    """
    ##--- calcuate histograms ---
    col_hist_dih, col_phi = d1.dihed_stat(sel1_kw, sel2_kw, sel3_kw, sel4_kw)
    ##--- save results ---
    np.savetxt(  fn_prefix+'dihed_CA5-C.dat', \
                 np.transpose( [col_phi, col_hist_dih]),\
                 fmt=['%f', '%f'], header='phi  P(phi)' )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


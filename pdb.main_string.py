from class_data import data
if __name__ == '__main__':
    import time as timer
    workdir = '/work/04201/tg835006/stampede2/LiuSum/Br/BrC5/430/'
    pdbfilename = 'vanhove430_2000.pdb' # vanhove_s calculation
    pdbfilename = workdir+pdbfilename
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = 'BrC5pdb_430_'
    fn_prefix = workdir+fn_prefix
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw) #, CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
    
    ##--- unwrap AN ---
    d1.unwrapall_AN()  # no need to unwrap if read from fix file and uxyz not deleted

    ##--- find fast AN ---
    fast_percentage = d1.find_AN_fast(16, 4.6, skip=1)  # interval_star= 16= 32000ps, r*=4.6A
    print('fast anion % = ' , fast_percentage)
    
    ##--- find string ---
    pns_histo = d1.find_AN_string(16, 2.5, maxlength=20, skip=1, include_rattle_ions = True)
    
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

    ##--- save
    np.savetxt(    fn_prefix+'pns_string.dat', np.transpose( np.vstack( (pns_histo[1], pns_histo[0]) )), fmt=['%d', '%f'], header='dist, P(ns)'      )


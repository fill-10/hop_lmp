from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = '../Tf2C2_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[sel.index%15==0]", 'type']
    CT_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[sel.index%15==2]", 'type']
    
    pdbfilename = '../../nvtprod_every10ps_00-10ns.pdb'
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    d1.allframes = d1.allframes[:-1]
    print('last timestep read: ', d1.allframes[-1].time)
    
    pdbfilename = '../../nvtprod_every10ps_10-20ns.pdb'
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    d1.allframes = d1.allframes[:-1]
    print('last timestep read: ', d1.allframes[-1].time)
    
    pdbfilename = '../../nvtprod_every10ps_30-40ns.pdb'
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    d1.allframes = d1.allframes[:-1]
    print('last timestep read: ', d1.allframes[-1].time)
    
    pdbfilename = '../../nvtprod_every10ps_40-50ns.pdb'
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    print('last timestep read: ', d1.allframes[-1].time)
    print('Number of frames:', len(d1.allframes))
    time_col, P1_col, P2_col = d1.RCF(10) #resolution = 10ps
    #--- save RCF 1 and 2 ---
    np.savetxt(    fn_prefix+'RCF_1.dat', np.transpose( np.vstack( (time_col, P1_col) )), fmt=['%f', '%f'], header='t/ps, C1'      )
    np.savetxt(    fn_prefix+'RCF_2.dat', np.transpose( np.vstack( (time_col, P2_col) )), fmt=['%f', '%f'], header='t/ps, C2'      )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


from class_data import data
if __name__ == '__main__':
    import time as timer
    #pdbfilename = '../nongauss_every100ps.pdb' # rdf calculation
    pdbfilename = '../nongauss_every10ps_0-10ns.pdb'
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = '../AmC2_'
    
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(11,11+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(1,1+10), 2, 40, 1, [30,1171], 15, "sel[ ((sel.index%30<=9 ) & (sel.index%30>=5)) | ( (sel.index%30>=15 ) & (sel.index%30<=24)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw )# , CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
 
    ##--- non gaussian parameter ---
    ##--- unwrap anion coordinates --- # may not necessary for some cases
    #d1.unwrapall_AN()  # do not unwrap if read from pdb. pdb has unwrapped data already.
    #print('unwrapped data')
    timestep_col, nongauss_col = d1.nongauss_AN_avg(1, 10, 1000) # start_interval = 1 frame,  100ps/frame, maxattemp =500 by default
    print('non gauss calculated')
    
    ##--- MSD

    timecol, msd_col  = d1.msd_AN_avg(1, 10, 1000)
    print('MSD calculated')


    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

    ##--- save non gaussian parameter
    np.savetxt(    fn_prefix+'nongauss_fine.dat', np.transpose( [timestep_col, nongauss_col] ), fmt=['%d', '%f'], header='time/ps, nongauss'      )
    np.savetxt(    fn_prefix+'MSD_fine.dat',      np.transpose( [timestep_col, msd_col]      ), fmt=['%d', '%f'], header='time/ps, msd'           )

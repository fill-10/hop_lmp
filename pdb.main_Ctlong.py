from class_data import data

if __name__ == '__main__':
    import time as timer
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = '../AmC2_'
 
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1 # save atoms and ions
    ##--- read ht from file ---
    d1.L_all_ht = np.loadtxt('./htlong.dat', dtype=int)
    print(d1.L_all_ht.shape[0], ' * ', d1.L_all_ht.shape[1])
    ##--- Ct ---
    time_col, Ct_col = d1.Ct( 100, 500 )
    np.savetxt(    fn_prefix+'Ctlong.dat', np.transpose( [time_col, Ct_col] ), fmt=['%d', '%f'], header='time/ps, Ct'      )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('Ct long calculated')
    print('time used in sec: ', stop-start)

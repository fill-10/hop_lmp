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
    d1.L_all_ht = np.loadtxt('./ht1.dat', dtype=int)
    d1.L_all_ht = np.append(d1.L_all_ht,  np.loadtxt('./ht2.dat', dtype=int)[1:], axis=0)
    d1.L_all_ht = np.append(d1.L_all_ht,  np.loadtxt('./ht3.dat', dtype=int)[1:], axis=0)
    d1.L_all_ht = np.append(d1.L_all_ht,  np.loadtxt('./ht4.dat', dtype=int)[1:], axis=0)
    print(d1.L_all_ht.shape[0], ' * ', d1.L_all_ht.shape[1])
    ##--- Ct ---
    time_col, Ct_col = d1.Ct( 10, 500 )
    np.savetxt(    fn_prefix+'Ct.dat', np.transpose( [time_col, Ct_col] ), fmt=['%d', '%f'], header='time/ps, Ct'      )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)

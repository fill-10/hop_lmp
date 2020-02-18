from class_data import data

if __name__ == '__main__':
    import time as timer
    ##--- output general settings ---=
    import numpy as np
    fn_prefix = '../Pf6C2_'
 
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=1 # save atoms and ions
    ##--- read ht from file ---
    d1.L_all_ht = np.loadtxt('./ht1.dat', dtype=int)
    d1.L_all_ht = np.append(d1.L_all_ht,  np.loadtxt('./ht4.dat', dtype=int)[1:], axis=0)
    d1.L_all_ht = np.append(d1.L_all_ht,  np.loadtxt('./ht3.dat', dtype=int)[1:], axis=0)
    d1.L_all_ht = np.append(d1.L_all_ht,  np.loadtxt('./ht4.dat', dtype=int)[1:], axis=0)
    print(d1.L_all_ht.shape[0], ' * ', d1.L_all_ht.shape[1])
    ##--- St ---
    time_col, St_col = d1.St( 10 )
    np.savetxt(    fn_prefix+'St_new.dat', np.transpose( [time_col, St_col] ), fmt=['%d', '%f'], header='time/ps, St'      )
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('St calculated.')
    print('time used in sec: ', stop-start)

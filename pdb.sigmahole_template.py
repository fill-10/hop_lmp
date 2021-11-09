from class_data import data
if __name__ == '__main__':
    import time as timer
    import numpy as np
    import pandas as pd 
    from pbc_dist import pbc_vec
    ##--- read-in the whole file ---
    #pdbfilename = '/Users/xuboluo/Documents/peptoid/luo_remake/Npe-4Br/antiparallel/test_trj.pdb'
    #pdbfilename = '/Users/xuboluo/Documents/peptoid/luo_remake/Npe-4Br/antiparallel_sigma/trj_000-004.pdb'
    pdbfilename = '/Users/xuboluo/Documents/peptoid/luo_remake/Npe-4Br/parallel/trj_000_224.pdb'
    pdbfilename = '/Users/xuboluo/Documents/peptoid/luo_remake/Npe-4Br/antiparallel/trj_000_399.pdb'
    pdbfilename = '/Users/xuboluo/Documents/peptoid/luo_remake/Npe-4Cl/parallel/trj_000_124.pdb'
    pdbfilename = '/Users/xuboluo/Documents/peptoid/luo_remake/Npe-4Cl/antiparallel/trj_000_399.pdb'
    ##--- output general settings ---
    #fn_prefix = '../Cl_para'
    fn_prefix = '../Cl_anti'
    ##--- timer start ---
    start = timer.perf_counter()
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0 #do not clear L_atoms

    ##--- read in pdb traj output ---
    # kw= [ molid, iontype, deg_poly, ion_per_mono, [iloc_terminus], Natom_ION, COM_col]
    # In this case, cannot use the ion_gen, please select atoms manually.
    d1.read_all_pdb(pdbfilename )
    print('last timestep read: ', d1.allframes[-1].time)
    for frame in d1.allframes:
        frame.time = frame.time*10
        frame.L_AN = frame.L_atom[frame.L_atom['type']=='CLH'].reset_index(drop=True)
        frame.L_CT = frame.L_atom[frame.L_atom['type']=='LP1'].reset_index(drop=True)
        del frame.L_atom
    d1.wrapall_L()
    ##--- calcuate histograms ---
    # find the coordinated LP1s for each BRH
    norm_hist_atom, norm_hist_mol = d1.find_asso_AN_CT(3.0)
    #print('last timestep read: ', d1.allframes[-1].L_AN, d1.allframes[-1].L_CT)
    # exclude the LP1 of the same monomer
    for frame in d1.allframes:
        L_asso_atom = [] # empty column
        for (idx_A, row_A) in frame.L_AN.iterrows():
            asso_atom = []
            for Ct in row_A['asso_atom']:
                if Ct == row_A['id']+1:
                    pass
                else:
                    asso_atom += [Ct]
            L_asso_atom += [asso_atom]
        frame.L_AN['asso_atom'] = L_asso_atom
    """
    # angle dist
    angle_hist = np.array([])
    # each snapshot
    for frame in d1.allframes:
        L_angle_at1 =   []
        L_angle_at2 =   []
        L_angle_at3 =   []
        for (idx, row) in frame.L_AN.iterrows():
            # build x y z matrix
            for i in range(0, len(row['asso_atom']) ):
                L_angle_at2 += [ [row['x'], row['y'], row['z']] ]  #center AN
                sg_self = frame.L_CT[ frame.L_CT['id'] == row['id']+1 ].iloc[0] #self CT as pd.Series
                L_angle_at1 += [ [ sg_self['x'], sg_self['y'], sg_self['z']] ]
                sg_out  = frame.L_CT[ frame.L_CT['id'] == row['asso_atom'][i]].iloc[0] # as pd.Series
                L_angle_at3 += [ [sg_out['x'],  sg_out['y'] , sg_out['z']  ] ]
        # vectors: LP-halogen-LP
        vec_self = np.array(L_angle_at1) - np.array(L_angle_at2)
        vec_out  = np.array(L_angle_at3) - np.array(L_angle_at2)
        # apply pbc
        pbc_vec(vec_out, [frame.deltaX, frame.deltaY, frame.deltaZ])
        pbc_vec(vec_self,[frame.deltaX, frame.deltaY, frame.deltaZ])
        # arccos
        cos_angle =  np.sum(np.multiply( vec_self, vec_out), axis=1 ) \
                  /  np.linalg.norm(vec_self, axis=1 ) \
                  /  np.linalg.norm(vec_out , axis=1 )
        L_angle = np.arccos( cos_angle)/np.pi*180
        c_hist, c_bins = np.histogram( L_angle, bins=np.arange(-0.5,180.5,1) )
        # add to all snapshots
        if angle_hist.shape[0]:
            angle_hist = np.vstack( (angle_hist, c_hist ) )
        else:
            angle_hist = c_hist
    sum_angle_hist = np.sum( angle_hist, axis =0 )
    avg_angle_hist = sum_angle_hist/np.sum(sum_angle_hist)
    ##--- save results ---
    np.savetxt(  fn_prefix+'hist_asso_atom.dat', \
                 np.transpose( [norm_hist_atom[1][:-1]-1, norm_hist_atom[0][:] ]),\
                 fmt=['%d', '%f'], header='n  P(n)' )
    np.savetxt(  fn_prefix+'hist_angle.dat', \
                 np.transpose( [np.arange(-0.5,180.5,1)[:-1]+0.5, avg_angle_hist] ),\
                 fmt=['%f', '%f'], header='alpha  P(alpha)' )
    """
    # inter connection:
    L_num_inter = []
    for frame in d1.allframes:
        interconnect =0
        L_asso_atom = [] # empty column
        for (idx_A, row_A) in frame.L_AN.iterrows():
                for  idx_B, row_B in frame.L_AN.iterrows():
                    if row_A['id']+1 in row_B['asso_atom'] and \
                        row_B['id']+1 in row_A['asso_atom']:
                            interconnect += 1
                            continue
        L_num_inter.append(interconnect/frame.L_AN.shape[0])
    print(np.mean(L_num_inter))
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


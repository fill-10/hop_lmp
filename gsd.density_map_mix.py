import numpy as np
import sys
import os
#parentdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..' ) )
#sys.path.append(parentdir)
from class_data import *
from class_oneframe import *

###--- Read traj
def prep_traj( fn_in,  L_Nchain=[], drop_rbody = True , pbc_whole= True  ):
    #L_Nchain = [[Nchain1, NAperCh1], [...], [...]]
    d0 = data(savemem=0)
    d0.read_all_gsd(fn_in, is_mass= True ) #print(d1.allframes[-1].L_atom)

    ###--- Data cleansing
    for frame in d0.allframes:
        # unit conversion nm to A
        frame.L_atom.loc[:, ['x', 'y', 'z' ] ] *= 10.
        frame.deltaX *= 10
        frame.deltaY *= 10
        frame.deltaZ *= 10
        frame.xhi *= 10
        frame.yhi *= 10
        frame.zhi *= 10
        frame.xlo *= 10
        frame.ylo *= 10
        frame.zlo *= 10
        #for vmd only, shift half box
        #frame.L_atom[ 'x' ] += frame.deltaX *0.5
        #frame.L_atom[ 'y' ] += frame.deltaY *0.5
        #frame.L_atom[ 'z' ] += frame.deltaZ *0.5

        # Repair Topology 
        # drop rigid body center of mass
        if drop_rbody:
            mask = ~frame.L_atom['type'].str.contains( 'rbody' )
            frame.L_atom = frame.L_atom[mask]
            frame.L_atom = frame.L_atom.reset_index( drop= True )
            frame.Natom = frame.L_atom.shape[0]
        # repair 'id' col
        frame.L_atom['id'] = frame.L_atom.index + 1
        # recover molid 
        col_mol = []
        acc_Nchain = 0
        Nmol = 0
        for t_chain in L_Nchain: # iterate each type of mol
            Nchain = t_chain[0]
            NAperCh = t_chain[1]
            Nmol += Nchain
            for i in range(acc_Nchain+0, acc_Nchain+Nchain):
                col_mol += [ i+1 for j in range(0,  NAperCh ) ]
            acc_Nchain  += Nchain
        frame.L_atom['mol' ] = col_mol

        # Make molecules whole and minimize ix iy iz
        if pbc_whole:
            for i in range(0, Nmol):
                molmask = frame.L_atom['mol'] == i + 1
                for col in ['ix', 'iy', 'iz']:
                    ix_res = frame.L_atom.loc[molmask, col].values.copy()
                    #idx_center =  np.argmin( np.abs(ix_res))
                    idx_center =  int ( ix_res.shape[0]/2 ) # middle atom
                    val_min_abs =  ix_res[ idx_center ] # + -
                    ix_res -= val_min_abs
                    frame.L_atom.loc[molmask, col] =  ix_res
    return d0


if __name__ == '__main__' :
    import time as timer
    start = timer.perf_counter()
    fn_in = '../50_50/rigid_slab.gsd'
    #fn_in = '../50_50/slab_init.gsd'
    # no need ot process pbc=_whole
    d_traj = prep_traj( fn_in, [[50, 148], [50, 196] ], drop_rbody = True , pbc_whole = False )
    # last 400 frames
    d_traj.allframes = d_traj.allframes[-400:]
    # remove COM translation
    d_traj.com_remove()
    # select tdp-43
    for frame in d_traj.allframes:
        mask = frame.L_atom['mol']<=50 # 1-50: tdp-43
        frame.L_AN = frame.L_atom.loc[ mask , : ]
    # density map, [-200, 200], bin size 1. Angstrom, dim=2: z direction
    res, bin_centers = d_traj.density_map_AN( -200. , 200., 1., dim = 2 )
    # mass unit: g/mol
    # density unit: g/mol / A^3
    # unit conversion: 1 mol: (6.02*10^23), 1 A^3 = 10^-24 cc
    res = res / 6.02 * 10  # g / cc
    # average +z and -z directions
    res = 0.5 * ( res + np.flip(res) )
    # save result
    np.savetxt( '../50_50/density_map_z_43.txt', np.vstack([bin_centers, res]).T, fmt=['%.2f', '%.6f'] )
    # select B8
    for frame in d_traj.allframes:
        mask = frame.L_atom['mol']>50 # 51-100: B8
        frame.L_AN = frame.L_atom.loc[ mask , : ]
    # density map, [-200, 200], bin size 1. Angstrom, dim=2: z direction
    res, bin_centers = d_traj.density_map_AN( -200. , 200., 1., dim = 2 )
    # mass unit: g/mol
    # density unit: g/mol / A^3
    # unit conversion: 1 mol: (6.02*10^23), 1 A^3 = 10^-24 cc
    res = res / 6.02 * 10  # g / cc
    # average +z and -z directions
    res = 0.5 * ( res + np.flip(res) )
    # save result
    np.savetxt( '../50_50/density_map_z_B8.txt', np.vstack([bin_centers, res]).T, fmt=['%.2f', '%.6f'] )

    stop  = timer.perf_counter()
    print('time used in sec: ', stop - start )


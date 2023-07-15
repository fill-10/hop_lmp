from class_data import data
from class_oneframe import oneframe

if __name__ == '__main__':
    import time as timer
    import numpy as np
    import pandas as pd 
    from pbc_dist import pbc_vec
    ##--- output general settings ---
    fn_prefix = '../'
    ##--- timer start ---
    start = timer.perf_counter()
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0 #do not clear L_atoms

    ##--- read in gro traj ---
    grofilename = '../test.gro'
    d1.read_all_gro(grofilename )
    d2 = data()
    d2.save_mem=0 #do not clear L_atoms
    grofilename = '../test_urea.gro'
    d2.read_all_gro(grofilename )
    for i in range( 0, len( d1.allframes ) ):
        d1.allframes[i].L_atom = \
            d1.allframes[i].L_atom.append( \
                d2.allframes[i].L_atom, ignore_index = True )
        d1.allframes[i].L_atom.reset_index(drop = True)
        d1.allframes[i].L_atom['id'] = d1.allframes[i].L_atom.index + 1
        d1.allframes[i].Natom = d1.allframes[i].L_atom.shape[0]
        #print(d1.allframes[i].L_atom)
    del d2 # save memory
    print('last timestep read: ', d1.allframes[-1].time)
    print( 'number of frames: ', len(d1.allframes) )
    ##--- move solvent ---
    Nchain = 12 * 4
    Napc = 685 # atoms per chain
    d1.prep_sasa( \
            "(self.L_atom['res'] == 'NDC') | (self.L_atom['res'] == 'NTE')",\
            " self.L_atom['res'] == 'UREA' ",                    0.3)


    """
    for i in range(0, len(d1.allframes)):
        print('frame: ', i)
        # set up cut-off for shielding:
        Rcut = 0.3 # 0.3 nm > 2*r_probe . 2*0.14
        d1.allframes[i].prep_sasa( \
            "(self.L_atom['res'] == 'NDC') | (self.L_atom['res'] == 'NTE')",\
            " self.L_atom['res'] == 'UREA' ",                    0.3        )

        # NOTE: Natom of each frame is different,
        # because of the selected urea for copy!
        # VMD cannot visualize the trajectory due to this difference,
        # but the tcl script can be executed with no problems.
        # Go!
    """
    """
        ##--- unit conversion ---
        frame.L_atom['xu'] *= 10
        frame.L_atom['yu'] *= 10
        frame.L_atom['zu'] *= 10
        frame.xlo *= 10
        frame.xhi *= 10
        frame.ylo *= 10
        frame.yhi *= 10
        frame.zlo *= 10
        frame.zhi *= 10
        frame.deltaX *= 10
        frame.deltaY *= 10
        frame.deltaZ *= 10
    """
    # export
    d1.export_all_gro(fn_prefix+'test_pep_urea.gro')
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


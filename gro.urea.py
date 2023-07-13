from class_data import data
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

    for i in range(0, len(d1.allframes)):
        # select Y direction:
        # be carefull! The selection is based on visuallization.
        sel_y_mol  = d1.allframes[i].L_atom.loc[ \
            (  (  d1.allframes[i].L_atom['yu']  > d1.allframes[i].deltaY *0.6 ) \
          &    (  d1.allframes[i].L_atom['res'] == 'UREA' )                     ) ,  'mol'].drop_duplicates( ).reset_index(drop=True)
        selrowy =  d1.allframes[i].L_atom['mol'].isin( sel_y_mol )
        # select Z direction
        sel_z_mol  = d1.allframes[i].L_atom.loc[ \
            (    (  d1.allframes[i].L_atom['zu']  < d1.allframes[i].deltaZ *0.4 ) \
              &  (  d1.allframes[i].L_atom['res'] == 'UREA' )                     ) , \
            'mol'                              ].drop_duplicates( ).reset_index(drop=True)
        selrowz =  d1.allframes[i].L_atom['mol'].isin( sel_z_mol )

        ##--- delete UREA and TFA to save time ---
        #    frame.L_atom = frame.L_atom[ \
        #    ( frame.L_atom['res'] == 'NDC' )  | \
        #    ( frame.L_atom['res'] == 'NTE' )  ]
        #
        ##-- move molecules
        # move Y
        d1.allframes[i].L_atom.loc[ selrowy, 'yu'] -= d1.allframes[i].deltaY
        # move Z
        d1.allframes[i].L_atom.loc[ selrowz, 'zu'] += d1.allframes[i].deltaZ

        ##-- copy X after moving moleculess for real UREA coverage by pbc
        # Do the selection after moving molecules in Y and Z.
        # 1/5 box should be enough, since there are 6 stacks, 4/5. 
        #
        # one side (left):
        sel_x_mol  = d1.allframes[i].L_atom.loc[ \
            (    (  d1.allframes[i].L_atom['xu']  >= d1.allframes[i].deltaX * 4 / 5 ) \
              &  (  d1.allframes[i].L_atom['res'] == 'UREA' )                     ) , \
            'mol'                              ].drop_duplicates( ).reset_index(drop=True)
        selrowx =  d1.allframes[i].L_atom['mol'].isin( sel_x_mol )
        # deep copy
        more_sol_l =   d1.allframes[i].L_atom.loc[ selrowx, : ].copy( deep = True )
        more_sol_l['xu'] -= d1.allframes[i].deltaX
        # the other side (right)
        sel_x_mol  = d1.allframes[i].L_atom.loc[ \
            (    (  d1.allframes[i].L_atom['xu']  <= d1.allframes[i].deltaX * 1 / 5 ) \
              &  (  d1.allframes[i].L_atom['res'] == 'UREA' )                     ) , \
            'mol'                              ].drop_duplicates( ).reset_index(drop=True)
        selrowx =  d1.allframes[i].L_atom['mol'].isin( sel_x_mol )
        #
        more_sol_r =   d1.allframes[i].L_atom.loc[ selrowx, : ].copy( deep = True )
        more_sol_r['xu'] += d1.allframes[i].deltaX

        # copy by append
        d1.allframes[i].L_atom = d1.allframes[i].L_atom.append(more_sol_l, ignore_index = True )
        d1.allframes[i].L_atom = d1.allframes[i].L_atom.append(more_sol_r, ignore_index = True )
        # reset id
        d1.allframes[i].L_atom.reset_index(drop = True)
        d1.allframes[i].L_atom['id'] = d1.allframes[i].L_atom.index + 1
        # update Natom
        d1.allframes[i].Natom = d1.allframes[i].L_atom.shape[0]
        # NOTE: Natom of each frame is different,
        # because of the selected urea for copy!
        # VMD cannot visualize the trajectory due to this difference,
        # but the tcl script can be executed with no problems.
        # Go!
        ##--- unit conversion ---
        """
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
        ##----------------------- 
    # export
    d1.export_all_gro(fn_prefix+'trj_P4_pep_urea.gro')
    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


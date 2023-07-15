from class_data import data

def pbcsolv(prot, sol, rcut, di, box_hi, box_lo, L_atom ):
# only for the solvent box of gmx trjconv -pbc whole,
# i.e., the box is [0, boxX(Y/Z)]
    # upper and lower bounds, 'di' is the direction,
    # which can be 'xu'/'yu'/'zu', etc. 
    # rcut is the cut off for the 'shield' beyond the bound.
    print(L_atom.shape)
    print('direction: ',  di )
    bup = prot[di].max() + rcut
    blo = prot[di].min() - rcut
    # empty DataFrame for copy:
    atom_copy = pd.DataFrame( columns=L_atom.columns )
    # extend above upper bound:
    if  bup > box_hi: # move or add sol mol to upper bound
        sel_mol = sol.loc[  sol[ di ] <=  bup - box_hi + box_lo ,  :  ]
        sel_molid = sel_mol['mol'].drop_duplicates( ).reset_index(drop=True)
        sel_row = L_atom['mol'].isin( sel_molid )
        # copy or move:
        print('avail space low side: ', blo-box_lo)
        print('req space hi side: ', bup - box_hi)
        if blo - box_lo > bup - box_hi : # True: move
            print(blo-box_lo)
            print('move to cover hi' )
            L_atom.loc[ sel_row, di ] += box_hi - box_lo
            pass
        else : # False: copy
            print('copy to cover hi' )
            more_atom = sol.loc[ sel_row, : ].copy( deep = True )
            more_atom.loc[:, di ]  +=  box_hi - box_lo
            #print( more_atom)
            atom_copy = atom_copy.append( more_atom, ignore_index = True )
            #print( atom_copy)
    # extend below lower bound:
    if blo < box_lo: #  to lower
        sel_mol = sol.loc[  sol[ di ] >=  blo + box_hi - box_lo ,  :  ]
        sel_molid = sel_mol['mol'].drop_duplicates( ).reset_index(drop=True)
        sel_row = L_atom['mol'].isin( sel_molid )
        # copy or move:
        print('avail space hi side: ', box_hi - bup)
        print('req space lo side: ', box_lo - blo)
        if box_hi - bup > box_lo - blo  : # True: move
            print('move to cover lo' )
            L_atom.loc[ sel_row, di ] -= box_hi - box_lo
            pass
        else : # False: copy
            print('copy to cover lo' )
            more_atom = sol.loc[ sel_row, : ].copy( deep = True )
            more_atom.loc[:, di ]  -=  box_hi - box_lo
            #print( more_atom)
            atom_copy = atom_copy.append( more_atom, ignore_index = True )
            #print( atom_copy)
    return  atom_copy

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
        print('frame: ', i)
        # set up cut-off for shielding:
        Rcut = 0.3 # 0.3 nm > 2*r_probe . 2*0.14
        # X
        dim = 0
        for direction in ['xu', 'yu', 'zu']:
            if dim == 0:
                box_lo = d1.allframes[i].xlo
                box_hi = d1.allframes[i].xhi
            elif dim == 1:
                box_lo = d1.allframes[i].ylo
                box_hi = d1.allframes[i].yhi
            elif dim == 2:
                box_lo = d1.allframes[i].zlo
                box_hi = d1.allframes[i].zhi
            print( box_lo, box_hi )
            more_copy = pbcsolv( \
                d1.allframes[i].L_atom.loc[ \
                    ( d1.allframes[i].L_atom['res'] == 'NDC' ) \
                  | ( d1.allframes[i].L_atom['res'] == 'NTE' ) , : ] ,  \
                d1.allframes[i].L_atom.loc[ \
                    d1.allframes[i].L_atom['res'] == 'UREA'    , : ] ,  \
                Rcut, direction, box_hi, box_lo , d1.allframes[i].L_atom )

            d1.allframes[i].L_atom = d1.allframes[i].L_atom.append( \
                more_copy, ignore_index = True )
            # counter +1 for x y z
            dim += 1
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


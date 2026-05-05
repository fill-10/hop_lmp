from class_data import *
def gsd2lmp_whole(fin, fout, topol_chain = [[100,148]], drop_rbody = True , pbc_whole= True ):
    # Load the GSD trajectory
    d = data()
    d.read_all_gsd(fin)
    for frame in d.allframes:
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
        frame.L_atom[ 'x' ] += frame.deltaX *0.5
        frame.L_atom[ 'y' ] += frame.deltaY *0.5
        frame.L_atom[ 'z' ] += frame.deltaZ *0.5

        # Repair Topology 
        # drop rigid body center of mass
        if drop_rbody:
            mask = ~frame.L_atom['type'].str.contains( 'rbody' )
            frame.L_atom = frame.L_atom[mask]
            frame.L_atom = frame.L_atom.reset_index( drop= True )
            frame.Natom = frame.L_atom.shape[0]
        # repair 'id' col
        frame.L_atom['id'] = frame.L_atom.index + 1
        # repair name, type
        frame.L_atom['res'] = frame.L_atom['type']
        frame.L_atom['type'] = ['C'] * frame.Natom
        # recover molid 
        col_mol = []
        Natom = 0
        Nmol = 0
        for mol in topol_chain:
            Nchain = mol[0]
            Nperchain = mol[1]
            # build new 'mol' col
            for i in range(0, Nchain):
                col_mol += [ i+1 + Nmol ] * Nperchain
            # counter
            Natom += Nchain * Nperchain
            Nmol += Nchain
        # check Natom dim
        if Natom != frame.L_atom.shape[0]:
            print('inconsistent [[Nchain, Nperchain]] with frame')
            return -1
        # assign 'mol' col
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

    d.export_all_lmptrj( fout )
    return 0

if __name__ == '__main__':
    # topol_chain param is consistent with drop_rbody, 
    # nojump
    """
    for i in range(0, 1):
        _ = gsd2lmp_whole(  f'../pure/rigid_slab{i}.gsd',
                            f'../pure/rigid_slab{i}.dump',
                            [  [ 100 , 148 ]  ],
                            drop_rbody = True ,
                            pbc_whole  = False )
    for i in range(0, 2):
        _ = gsd2lmp_whole(  f'../pure_v2/rigid_slab{i}.gsd',
                            f'../pure_v2/rigid_slab{i}.dump',
                            [  [ 100 , 148 ]  ],
                            drop_rbody = True ,
                            pbc_whole  = False )
    for i in range(0, 1):
        _ = gsd2lmp_whole(  f'../pure/rigid_slab{i}.gsd',
                            f'../pure/rigid_slab{i}_whole.dump',
                            [  [ 100 , 148 ]  ],
                            drop_rbody = True ,
                            pbc_whole  = True )

    for i in range(0, 2):
        _ = gsd2lmp_whole(  f'../pure_v2/rigid_slab{i}.gsd',
                            f'../pure_v2/rigid_slab{i}_whole.dump',
                            [  [ 100 , 148 ]  ],
                            drop_rbody = True ,
                            pbc_whole  = True )
    """
    _ = gsd2lmp_whole(  f'../50_50_v2/rigid_slab.gsd',
                        f'../50_50_v2/rigid_slab_whole.dump',
                        [   [ 50 , 148 ],
                            [ 50 , 196 ]  ],
                        drop_rbody = True ,
                        pbc_whole  = True )
    """
    _ = gsd2lmp_whole(  f'../75_25_v2/rigid_cube_em.gsd',
                        f'../75_25_v2/rigid_cube_em.dump',
                        [   [ 75 , 148 ],
                            [ 25 , 196 ]  ],
                        drop_rbody = True ,
                        pbc_whole  = False )

    """

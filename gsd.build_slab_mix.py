from class_data import *
import gsd 
# Load the GSD trajectory
def make_slab(fin_cube, fout_slab, topol_chain = [[100,148]], ext_dim = ['z'], ext_ratio=[10.]):
    d = data()
    #d.read_all_gsd('../1rigid/rigid_cube_final.gsd')
    d.read_all_gsd(fin_cube)
    frame = d.allframes[-1]
    frame.L_atom['id'] = frame.L_atom.index
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
    frame.L_atom['mol' ] = col_mol # new 'mol' col
    print(frame.L_atom)
    if Natom != frame.L_atom.shape[0]:
        print('inconsistent [[Nchain, Nperchain]] with frame')
        return -1
    # make slab
    # make molecules whole and wrap
    # x & y: minimize abs(ix) abs(iy)
    # z: minimize abs(iz) and unwrap
    newx = []
    newy = []
    newz = []
    for i in range(0, Nmol) :
        print(i)
        molmask = frame.L_atom['mol'] == i + 1
        for col in ['ix', 'iy', 'iz']:
            ix_res = frame.L_atom.loc[molmask, col].values.copy()
            #idx_center =  np.argmin( np.abs(ix_res))
            idx_center =  int ( ix_res.shape[0]/2 ) # middle atom
            val_min_abs =  ix_res[ idx_center ] # + -
            ix_res -= val_min_abs
            frame.L_atom.loc[molmask, col] =  ix_res
        for col in ext_dim :
            ix_res = frame.L_atom.loc[molmask, 'i'+ col].values.copy()
            box_dim = eval (' frame.delta' + col.upper() )
            exec ( "new" + col + ".append(  frame.L_atom.loc[ molmask, col ] + ix_res * box_dim )" )
    # save new z
    for col in ext_dim :
        frame.L_atom[ col ] = np.hstack( eval ( 'new'+col) ).reshape(-1 )
        frame.L_atom[ 'i'+col ] = 0
    # new box dim  for slab
    for i, col in enumerate(ext_dim):
        exec ( "frame."+col+"hi = " + "frame."+col+"hi * " + str (ext_ratio[i] ) )
        exec ( "frame."+col+"lo = " + "frame."+col+"lo * " + str (ext_ratio[i] ) )
        exec ( "frame.delta"+col.upper()+"*= " + str (ext_ratio[i] ) )

    # save lammps trj for visual inspection
    d.export_all_lmptrj( '../50_50/slab_init.dump' )

    # save new slab gsd 
    traj = gsd.hoomd.open(fin_cube, 'r')
    sn = traj[-1]
    sn.particles.position = frame.L_atom.loc[:, ['x', 'y', 'z' ] ].values
    sn.particles.image = frame.L_atom.loc[:, ['ix', 'iy', 'iz' ] ].values
    sn.configuration.box = np.array( [frame.deltaX, frame.deltaY, frame.deltaZ , 0., 0., 0. ] ) 
    sn.configuration.step = 0 # reset timestep
    sn.validate() # validate xyz
    with gsd.hoomd.open(name = fout_slab, mode = 'w' ) as f:
        f.append(sn)
    # confirm info in the new gsd
    traj = gsd.hoomd.open(fout_slab, 'r')
    sn = traj[-1]
    print(sn.particles.position)
    print(sn.particles.image)
    print(sn.particles.velocity)
    print(sn.configuration.box)
    print(sn.configuration.step)
    return 0

if __name__ == '__main__' :
    #for cycle in range(0, 8):
    #    _ = make_slab(  '../pure/rigid_cube'+str(cycle)+'.gsd',
    #                    '../pure/slab_init'+str(cycle)+'.gsd' ,
    #                    50 , [ 'z' ] , [ 10. ]                      )
    _ = make_slab(  '../50_50_v2/rigid_cube_npt_final.gsd',
                    '../50_50_v2/slab_init.gsd' ,
                    [ [50, 148+1], [50, 196+2  ]] , [ 'z' ] , [ 10. ]    )

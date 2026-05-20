import numpy as np
import sys
import os
#parentdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..' ) )
#sys.path.append(parentdir)
from class_data import *
from class_oneframe import *
from pbc_dist import pbc_vec, pbc_x
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.sparse.csgraph import connected_components

###--- Read traj
def prep_traj( fn_in,  L_Nchain=[], N_last = 110, N_every = 1 ):
    #L_Nchain = [[Nchain1, NAperCh1], [...], [...]]
    d0 = data(savemem=0)
    d0.read_all_gsd(fn_in) #print(d1.allframes[-1].L_atom)
    d0.allframes = d0.allframes[ -N_last : : N_every] # discard first 1 us
    Nframe = len ( d0.allframes )
    print("number of frames: ", Nframe )
    #print(d2.allframes[-1].L_atom.shape)
    # unwrap
    d0.unwrapall_L()
    
    ###--- Data cleansing
    for frame in d0.allframes:
        #frame.L_atom['res'] = frame.L_atom['type']
        #frame.L_atom['type'] = ['C'] * frame.Natom
        # nm to A
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

        # pbc whole
        for i in range(0, Nmol):
            molmask = frame.L_atom['mol'] == i + 1
            for col in ['ix', 'iy', 'iz']:
                ix_res = frame.L_atom.loc[molmask, col].values.copy()
                #idx_center =  np.argmin( np.abs(ix_res))
                idx_center =  int ( ix_res.shape[0]/2 ) # middle atom
                val_min_abs =  ix_res[ idx_center ] # + -
                ix_res -= val_min_abs
                frame.L_atom.loc[molmask, col] =  ix_res

        # recover atom id
        frame.L_atom['id'] = frame.L_atom.index + 1
        # Do not drop rigid body center of mass
        #mask = ~frame.L_atom['type'].str.contains( 'rbody' )
        
        # pbc whole
        for i in range(0, Nmol):
            molmask = frame.L_atom['mol'] == i + 1

    return d0

def find_stack( d1, fn_out, L_Nchain=[] ):
    for frame in d1.allframes:
        # collect COM for each stack
        com_mask = frame.L_atom['type'].str.contains( 'rbody' )
        frame.L_AN = frame.L_atom[ com_mask ].reset_index(drop = True )
        
        # pca0 and pca2 for each rbody
        idx_body0 = \
               [ i for i in range(283-267,287-267+1)] \
            +  [ i for i in range(294-267,299-267+1)] \
            +  [ i for i in range(302-267,304-267+1)] \
            +  [ i for i in range(311-267,319-267+1)] \
            +  [ i for i in range(325-267,329-267+1)] \
            +  [ i for i in range(333-267,343-267+1)] \
            +  [ i for i in range(346-267,348-267+1)] \
            +  [ i for i in range(354-267,358-267+1)]
        idx_body0 = (np.array(idx_body0) + 1 ) #.tolist()

        # find local xyz for consitituent particles in each rbody
        acc_Nchain = 0
        L_pc0 = []
        L_pc2 = []
        for t_chain in L_Nchain: # iterate each type of mol
            Nchain = t_chain[0]
            NAperCh = t_chain[1]
            for i in range(acc_Nchain+0, acc_Nchain+Nchain):
                rbody_mask = i * NAperCh + idx_body0
                sel = frame.L_atom.iloc [ rbody_mask , : ]
                com = frame.L_AN.iloc[ [i] , :]
                sel = sel[ ['xu', 'yu', 'zu' ]].values
                com = com[ ['xu', 'yu', 'zu' ]].values
                # subtract com
                sel = sel - com
                # check if out of box
                if np.any( np.abs(sel) - np.array( [frame.deltaX/2, frame.deltaY/2, frame.deltaZ/2] ) > 0 ):
                    print('outofbox')
                    return -1 # find pc0 and pc2 , pc1 is useless pca = PCA(n_components = 3 )
                pca = PCA( n_components = 3 )
                pca.fit( sel )
                pc0 = pca.components_[0] # the res pc0 is normalized
                pc2 = pca.components_[2]
                L_pc0.append( pc0 )
                L_pc2.append( pc2 )
            acc_Nchain  += Nchain
        # save all pc0 and pc2 xyz in L_AN
        frame.L_AN.loc[ :, ['pc0x', 'pc0y', 'pc0z'] ] = L_pc0
        frame.L_AN.loc[ :, ['pc2x', 'pc2y', 'pc2z'] ] = L_pc2
        frame.L_AN = frame.L_AN[['id', 'type', 'mol', 'x', 'y', 'z', 'ix', 'iy', 'iz', 'pc0x', 'pc0y', 'pc0z', 'pc2x', 'pc2y', 'pc2z']]
        # delete useless L_atom
        frame.L_atom = pd.DataFrame(columns=['id', 'type', 'mol', 'x', 'y', 'z', 'ix', 'iy', 'iz'] )
    
    ###--- pairwise contact map for rbody ---###
    Nrbody = frame.L_AN.shape[0]
    #print( Nrbody)
    for frame in d1.allframes:
        dist_mat = np.zeros( [Nrbody, Nrbody] , dtype = float )
        pc0_mat = np.zeros( [Nrbody, Nrbody] , dtype = float )
        pc2_mat = np.zeros( [Nrbody, Nrbody] , dtype = float )
        for idx , row in frame.L_AN.iterrows():
            # dist criterion
            ref_com = row[['x', 'y', 'z']].values
            L_vec = frame.L_AN[ ['x', 'y', 'z'] ] - ref_com 
            # pbc dist for com
            L_vec = pbc_vec( L_vec, [frame.deltaX, frame.deltaY, frame.deltaZ] ) 
            L_norm = np.linalg.norm( L_vec.astype(np.float32), axis = 1 ) # dist
            dist_mat[idx, :] = L_norm
            #
            # angle criterion
            ref_pc0 = row[['pc0x', 'pc0y', 'pc0z']].values
            ref_pc2 = row[['pc2x', 'pc2y', 'pc2z']].values
            pc0_mat[idx , : ]  = np.abs( np.dot( frame.L_AN[ ['pc0x', 'pc0y', 'pc0z'] ] , ref_pc0 ) )
            pc2_mat[idx , : ]  = np.abs( np.dot( frame.L_AN[ ['pc2x', 'pc2y', 'pc2z'] ] , ref_pc2 ) )
        #print(np.min(dist_mat[dist_mat>0.01]) )
        #print(pc0_mat)
        cmap_bin =  (
            ( dist_mat < 10.0 )   &
            ( pc0_mat > 0.9396 ) & 
            ( pc2_mat > 0.9659 )
                ).astype(int)  # cos(20) = 0.9396 cos(15)=0.9659
        #print(cmap_bin)

        num_clusters, labels = connected_components (
            csgraph = cmap_bin  ,
            directed = False ,
            return_labels = True , 
                                                    )
        frame.L_AN['mol'] = labels
        #print(frame.L_AN['mol'] )
    # export
    d1.export_all_lmptrj(fn_out)
    return 0

    #save 
    #np.savetxt( fn_out_pre + '_H.txt', H, fmt='%.6f' )

if __name__ == '__main__' :
    import time as timer
    start = timer.perf_counter()
    #d1.read_all_lmp( 'rigid_slab.dump')
    #L_fn_in = [ f'../2rigid_beta_only_3pb/rigid_slab{i}_uconn.gsd' for i in range(0,8) ]
    fn_in = '../pure_v2/rigid_slab1.gsd'
    fn_out = '../pure_v2/trj_slab1_stack.dump'
    d_traj = prep_traj( fn_in, [[100, 148+1]], 1300 , 1)
    _ = find_stack(d_traj, fn_out,  [[100, 148+1]] )
    stop  = timer.perf_counter()
    print('time used in sec: ', stop - start )

    """
    # Gemini vectorization ...
    # Using broadcasting: shape (N, 1, 3) - shape (1, N, 3) -> (N, N, 3)
    # Taking norm along axis=2 results in a pairwise (N, N) distance matrix
    diff_centroids = centroids[:, np.newaxis, :] - centroids[np.newaxis, :, :]
    dist_matrix = np.linalg.norm(diff_centroids, axis=2)
    
    # 2. Vectorized Absolute Dot Products for Orientation Alignment
    # Matmul of shape (N, 3) x (3, N) yields an (N, N) matrix containing all dot products
    pc1_dot_matrix = np.abs(np.dot(pc1_axes, pc1_axes.T))
    pc3_dot_matrix = np.abs(np.dot(pc3_axes, pc3_axes.T))
    """

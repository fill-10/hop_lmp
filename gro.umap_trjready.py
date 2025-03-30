from class_data import data
from class_oneframe import oneframe
import matplotlib.pyplot as plt
import umap
from timer import my_timer
import numpy as np
import pandas as pd 
import os
import pickle

if __name__ == '__main__':
    @my_timer
    def run(): # train umap model with 4M thf data
        ##--- initialize a data object ---
        d1 = data()
        d1.save_mem=0 #do not clear L_atoms

        ##--- read in gro traj ---
        # 3001 frames of rest2 isolated
        # 202  frames of dimer stack, c side
        # 202  frames of dimer stack, o side
        # 306  frames of 6layer c side
        # 144 frames of nanosheet
        # [0,3001) rest2, [3001,3203) stack2 c,  [3203, 3405) stack2 o, [3405, 3711) 6layers c, [3711, 3855) sheet
        grofilename = 'trj_rest2_2_6layers_inf_NDC_noh.gro'
        d1.read_all_gro(grofilename )
        # confirm total number of frames
        print ( 'number of all  frames: ', len(d1.allframes) )
        # trim to Ndc bb only
        for frame in d1.allframes:
            #sel_row = ( frame.L_atom['res'] == 'NDC' ) & ( frame.L_atom['type'].str.get(0) != 'H' )
            sel_row = ( frame.L_atom['res'] == 'NDC' ) \
                &     (   ( frame.L_atom['type']== 'N') | ( frame.L_atom['type']== 'C' ) \
                        | ( frame.L_atom['type']== 'O') | ( frame.L_atom['type']== 'CA' ) \
            #            | ( frame.L_atom['type']== 'CD') | ( frame.L_atom['type']== 'CT' ) \
                      )
            #sel_row = ( frame.L_atom['res'] == 'NDC' ) \
            #    &     ( frame.L_atom['type']== 'N'   )
            frame.L_atom = frame.L_atom[sel_row].reset_index( drop = True )
            frame.Natom = frame.L_atom.shape[0]
        # gen zmat
        d1.gen_zmat(is_sqrt= True )
        # store zmat into one array
        myframes = np.array( [d1.allframes[0].L_zmat ])
        print ( 'shape of zmat:',  myframes.shape )
        for i in range (1 , len(d1.allframes) ):
            myframes = np.concatenate( [ myframes, [d1.allframes[i].L_zmat] ] ,  axis = 0 )
        #rescale with mean and std
        mymean = np.mean(myframes, axis = 0 )
        mystd = np.std(myframes, axis = 0 )
        myframes = ( myframes - mymean ) / mystd
        print (myframes.shape)
        # initialize umap obj: 'reducer' with default hypoparameters
        reducer = umap.UMAP(n_neighbors= 15 , min_dist=0.1 , random_state = 42 )
        # train umap, save transformed data in 'thf4m_train'
        thf4m_train = reducer.fit_transform( myframes  )
        # open a new dir for the hyperparam set
        new_dir = 'zmat_umap_NDCbb_15_0d1/'
        try :
            os.mkdir(new_dir)
        except :
            pass
        # save trained obj as a pickle file
        with open(new_dir+'umap_model.pkl', 'wb') as file:
            pickle.dump( reducer, file )
        # save tranformed data for 4M thf
        np.savetxt( new_dir+'umap_train.txt', thf4m_train , fmt = ['%.4f', '%.4f'] )
        # plot 4M thf umap
        x_re = thf4m_train[:3001, 0]
        y_re = thf4m_train[:3001, 1]
        x_2c = thf4m_train[3001:3001+202, 0]
        y_2c = thf4m_train[3001:3001+202, 1]
        x_2o = thf4m_train[3203:3203+202, 0]
        y_2o = thf4m_train[3203:3203+202, 1]
        x_6c = thf4m_train[3405:3405+306, 0]
        y_6c = thf4m_train[3405:3405+306, 1]
        x_sh = thf4m_train[3711:, 0]
        y_sh = thf4m_train[3711:, 1]
        plt.scatter( x_2c, y_2c, facecolors='none', edgecolors='green' , linewidth = .5 , s = 20, label = '2S C'  ) # 2S c
        plt.scatter( x_2o, y_2o, facecolors='none', edgecolors='orange',linewidth = .5 , s = 20, label = '2S O' ) # 2S o
        plt.scatter( x_6c, y_6c, facecolors='none', edgecolors='red', linewidth = .5 , s = 20, label = '6L C') # 6L
        plt.scatter( x_sh, y_sh, facecolors='none', edgecolors='black', linewidth = .5 , s = 20, label = 'sheet') # 6L
        #
        plt.scatter( x_re, y_re, facecolors='none', edgecolors='cyan', marker='D', s = 10, label = '4M thf' ) # sheet

        """
        #load pure thf rest
        del d2
        del d2zmat
        d2 = data()
        d2.save_mem=0 #do not clear L_atoms
        ##--- read in gro traj ---
        grofilename = '/Users/xuboluo/Documents/test_charmm/Ndc10-Nte10/Ac_thf_ld/rest/T300/trj_P2.gro'
        d2.read_all_gro(grofilename )
        d2.allframes=d2.allframes[::10]

        d3 = data()
        d3.save_mem=0 #do not clear L_atoms
        ##--- read in gro traj ---
        grofilename = '/Users/xuboluo/Documents/test_charmm/Ndc10-Nte10/Ac_thf_ld/rest/T300/trj_P3.gro'
        d3.read_all_gro(grofilename )
        d3.allframes=d3.allframes[::10]
        d2.allframes += d3.allframes[1:]
        del d3

        d3 = data()
        d3.save_mem=0 #do not clear L_atoms
        ##--- read in gro traj ---
        grofilename = '/Users/xuboluo/Documents/test_charmm/Ndc10-Nte10/Ac_thf_ld/rest/T300/trj_P4.gro'
        d3.read_all_gro(grofilename )
        d3.allframes=d3.allframes[::10]
        d2.allframes += d3.allframes[1:]
        del d3
        # confirm total number of frames
        print ( 'number of all  frames: ', len(d2.allframes) )

        # trim to Ndc bb only
        for frame in d2.allframes:
            #sel_row = ( frame.L_atom['res'] == 'NDC' ) & ( frame.L_atom['type'].str.get(0) != 'H' )
            sel_row = ( frame.L_atom['res'] == 'NDC' ) \
                &     (   ( frame.L_atom['type']== 'N') | ( frame.L_atom['type']== 'C' ) \
                        | ( frame.L_atom['type']== 'O') | ( frame.L_atom['type']== 'CA' ) \
            #            | ( frame.L_atom['type']== 'CD') | ( frame.L_atom['type']== 'CT' ) \
                      )
            #sel_row = ( frame.L_atom['res'] == 'NDC' ) \
            #    &     ( frame.L_atom['type']== 'N'   )
            frame.L_atom = frame.L_atom[sel_row].reset_index( drop = True )
            frame.Natom = frame.L_atom.shape[0]
        # zmat
        d2.gen_zmat(is_sqrt= True )
        # store zmat into one array
        d2zmat = np.array( [d2.allframes[0].L_zmat ])
        for i in range (1 , len(d2.allframes) ):
            d2zmat = np.concatenate( [ d2zmat, [d2.allframes[i].L_zmat] ] ,  axis = 0 )
        #rescale
        d2mean = np.mean(d2zmat, axis = 0 )
        d2std = np.std(d2zmat, axis = 0 )
        d2zmat = ( d2zmat - mymean ) / mystd
        print ( 'shape of d2.zmat:',  d2zmat.shape )
        # umap
        thf_test = reducer.transform( d2zmat )
        # plot
        x2 = thf_test[:, 0]
        y2 = thf_test[:, 1]

        plt.scatter( x2, y2, facecolors='none', edgecolors='olive', s = 10) # sheet
        
        np.savetxt( new_dir+'umap_thf.txt', thf_test , fmt = ['%.4f', '%.4f'] )
        """
        plt.legend(framealpha=0.)
        plt.xlabel('umap 1')
        plt.ylabel('umap 2')
        plt.show()



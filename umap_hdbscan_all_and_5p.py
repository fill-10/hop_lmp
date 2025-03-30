import numpy as np
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.cluster import HDBSCAN
from class_data import data
from class_oneframe import oneframe
import ase
import ase.visualize


def add_scatter_cluster_2d(x, y, labels, color_tab = mcolors.TABLEAU_COLORS):
    label_min = min(labels)
    label_max = max(labels)
    print( 'all labels: ', list(set(labels.astype(int) )) )
    colors = list(color_tab)
    Ncolors = len(colors)
    markers = ['D', 'h', '*', 'X', 's', 'v', '^']
    
    for i in range( 0, label_max+1):
        # indices of value i in label
        idx = [ j for j, k in enumerate(labels) if k == i]
        # x and y for cluster i
        x_i = x[idx]
        y_i = y[idx]
        # set color, i starts from 1, colors starts from 0
        color = colors[ (i) % Ncolors ]
        symbol = markers[ i // Ncolors ]
        #print('cluster: ', i,'color: ',  color)
        # add scatter plot
        plt.scatter( x_i, y_i, facecolors='none', edgecolors=color , marker = symbol , linewidth = .5, s = 20, label= str(i) )

def vis_avg_clsuter_2d(i, labels, data ): # avg of ith cluster
    # indices of value i in label
    indices = [ j for j, k in enumerate(labels) if k == i ]
    # select the frames for cluster i
    fms_sel = [  data.allframes[l] for l in indices ]
    print('num of cluster '+ str(i) + ' :', len(fms_sel) )
    # avg for xyz after zmat rebuild
    frame_avg = fms_sel[0]
    for m in range( 1 , len(fms_sel) ):
        frame = fms_sel[m]
        frame_avg.L_atom['xu'] = frame_avg.L_atom['xu'] + frame.L_atom['xu']
        frame_avg.L_atom['yu'] = frame_avg.L_atom['yu'] + frame.L_atom['yu']
        frame_avg.L_atom['zu'] = frame_avg.L_atom['zu'] + frame.L_atom['zu']
    # divide
    frame_avg.L_atom[ ['xu','yu','zu'] ] /= len(fms_sel)
    # visualize averge for cluster i using ase
    trj = []
    trj.append ( \
        ase.Atoms( frame_avg.L_atom['type'].str.get(0).to_list(),\
                   frame_avg.L_atom[['xu', 'yu', 'zu']].values *10 ) \
               ) # nm *10 -> angstrom
    ase.visualize.view( trj )

##-- set up data folder --##
sel='bb'
#umap_set='15_0d1'
umap_set='15_0d0'
#folder = 'zmat_umap_NDC'+sel+'_'+umap_set + '/'
folder = 'train_all_and_5p_umap_NDC'+sel+'_'+umap_set + '/'

##-- hdbscan for 4M thf --##
species = [ '4Mthf', 'water', '5pthf', 'purethf' ]
for s in species:
    fn = folder + 'umap_'+s+'.txt'
    reduced = np.loadtxt( fn )
    x = reduced[:, 0]
    y = reduced[:, 1]
    x_rest=x[     : 3001]
    y_rest=y[     : 3001]

    if s != 'purethf' and s != '5pthf' and s !='5p' :
        x_2c = x[3001 : 3203]
        y_2c = y[3001 : 3203]
        x_2o = x[3203 : 3405]
        y_2o = y[3203 : 3405]
        x_6c = x[3405 : 3711]
        y_6c = y[3405 : 3711]
        x_sh = x[3711 : ]
        y_sh = y[3711 : ]
    else:
        x_2c =  [ ]
        y_2c =  [ ]
        x_2o =  [ ]
        y_2o =  [ ]
        x_6c = x[3001 : 3307]
        y_6c = y[3001 : 3307]
        x_sh = x[3307 : ]
        y_sh = y[3307 : ]

    ##--- hdb scan ---
    hdb = HDBSCAN( min_cluster_size = 10 , min_samples = 3 , cluster_selection_epsilon = 0.2 )
    hdb.fit(reduced[:3001,:])
    print ( hdb.labels_, hdb.labels_.shape )
    ##--- set up font and size of plot ---##
    plt.rcParams['font.size'] = 25
    plt.rcParams['font.family'] = 'sans'
    plt.figure( figsize = (8,6) )
    # plot circles for dimer 6c and sheet
    plt.scatter( x_2c, y_2c , facecolors='none', edgecolors='green'  , marker = 'o', linewidth = .5, s = 29, )
    plt.scatter( x_2o, y_2o , facecolors='none', edgecolors='orange' , marker = 'o', linewidth = .5, s = 29, )
    plt.scatter( x_6c, y_6c , facecolors='none', edgecolors='red'    , marker = 'o', linewidth = .5, s = 29, )
    plt.scatter( x_sh, y_sh , facecolors='none', edgecolors='black'  , marker = 'o', linewidth = .5, s = 29, )
    # call function to plot all clusters
    add_scatter_cluster_2d(x, y, hdb.labels_)
    # plot axes
    plt.xlabel('umap 1')
    plt.ylabel('umap 2')
    #plt.axis( [ 7.5,17.5, -7,12])
    plt.tick_params( axis='both', direction='in', length=6)
    plt.tight_layout()
    #plt.legend(framealpha=0.)
    # show/save
    plt.show()
    #plt.savefig( folder + 'cluster_map_'+ s + '.png', dpi=600, transparent=True)
    #plt.cla()


    #np.savetxt( folder+'cluster_4Mthf.txt',\
    #            np.concatenate( [reduced, np.array( [ hdb.labels_] ).T], axis = 1 ), \
    #            fmt = ['%.4f', '%.4f', '%d'] \
    #        )

    ##-- visualiz avg cluster
    d1=data()
    d1.save_mem=0 #do not clear L_atoms
    if s == '4Mthf':
        grofilename = 'trj_rest2_2_6layers_inf_NDC_noh.gro'
    elif s == 'water' :
        grofilename = 'trj_water_rest_2_6_inf_NDC_noh.gro'
    elif s == '5pthf' or s =='5p' :
        grofilename = 'trj_5pthf_rest_6_inf_NDC_noh.gro'
    elif s == 'purethf' :
        grofilename = 'trj_thf_rest_6_inf_NDC_noh.gro'
    else:
        continue
    print('system: ', s )
    d1.read_all_gro(grofilename )
    #
    d1.allframes = d1.allframes[:3001]
    # trim to Ndc bb
    for frame in d1.allframes:
        sel_row = ( frame.L_atom['res'] == 'NDC' ) \
            &     (   ( frame.L_atom['type']== 'N') | ( frame.L_atom['type']== 'C' ) \
                    | ( frame.L_atom['type']== 'O') | ( frame.L_atom['type']== 'CA') )
        frame.L_atom = frame.L_atom[sel_row].reset_index( drop = True )
        frame.Natom = frame.L_atom.shape[0]
    # gen zmat
    d1.gen_zmat(is_sqrt= True )
    # rebuild xu yu zu to remove rot and trans
    d1.zmat2xyz()
    # visualize each cluster using input
    counter = 10
    while counter>=0:
        counter -= 1
        cluster_id = input('enter cluster id:')
        cluster_id = int(cluster_id)
        print('ith cluster: ', cluster_id)
        if cluster_id <0:
            break
        vis_avg_clsuter_2d(cluster_id, hdb.labels_[:3001], d1) # rest only


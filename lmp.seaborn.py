import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
#parentdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..' ) )
#sys.path.append(parentdir)
from class_data import *
from class_oneframe import *

def violin_df_prep(trj_obj):

    #Nframe = len(trj_obj.allframes)
    
    # find global max of stack size
    global_max_size = 1
    for frame in trj_obj.allframes:
        labels = frame.L_atom['mol']
        _, sizes = np.unique( labels, return_counts = True )
        global_max_size = max (global_max_size , np.max(sizes, axis = None ) )
    
    # histogram for each frame
    bin_edges = np.arange(0.5 , global_max_size + 1.5 , 1.0 )
    
    possible_sizes = np.arange(1, global_max_size + 1 )
    
    L_hist = []
    L_stat = []
    L_reconstr = []
    for frame_idx, frame in enumerate(trj_obj.allframes):
        labels = frame.L_atom['mol']
        _, sizes = np.unique( labels, return_counts = True )
        # density hist
        hist_density , _ = np.histogram( sizes, bins = bin_edges, density = True )
        # avg
        avg_size = np.mean(sizes)
        # std
        std_size = np.std(sizes)
        # save in list
        L_hist.append( hist_density.astype(np.float32) )
        L_stat.append( [ avg_size, std_size ] )

        # reconstructed data frame for plot
        for size in sizes:
            L_reconstr.append( [ frame_idx, size ] )
    
    L_hist = np.array(L_hist)
    L_stat = np.array(L_stat)
    df = pd.DataFrame( L_reconstr, columns = ['t', 'N'] )

    return df, L_hist, L_stat

def violin_plot(df,  fout_png, violin_width= 0.5 , sys_name = 'core' ):
    # plot
    fig, ax = plt.subplots( figsize = (16, 12))
    sns.set_theme(style='whitegrid' )

    # grid line
    #ax.grid( axis = 'y', linestyle='--', alpha=0.3, color='black', zorder = 0 )

    # Gnerate the timeline violins

    sns.violinplot( 
        data = df,
        x = df.columns[0], 
        y = df.columns[1],
        color = 'orange',
        inner = 'box', 
        width = violin_width,
        bw_adjust = 0.5 ,
        cut = 0,
        ax = ax 
    )
    
    # set axis 
    min_y = 0.8
    max_y = int ( df['N'].max() )
    ax.set_yticks(np.arange(0, max_y + 1, 5) )
    ax.set_xlabel('t', fontsize=48, labelpad=10)   
    ax.set_ylabel('N', fontsize=48, labelpad=10)   
    ax.tick_params(axis = 'both', which = 'major', labelsize = 36)


    # transparent background
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)

    plt.tight_layout()

    if fout_png:
        plt.savefig(fout_png, dpi=450, transparent= True )
    else :
        #plt.show()
        pass
    plt.close()

    # return df for future avg
    df.columns = ['sys', 'N']
    df['sys'] = sys_name
    return df

if __name__ == '__main__':
    import time as timer
    start = timer.perf_counter()
    # 7PY2
    d1=data(savemem = 0 )
    d1.read_all_lmp( '../pure_v2/trj_slab0_stack.dump' )
    d1.allframes = d1.allframes[ -200 :  : 2 ]
    df, _, stats =  violin_df_prep(d1 )
    #print(stats)
    df_7PY2 = violin_plot(df, None , sys_name='7PY2') #'seaborn_violin0.png'

    # 8CG3
    d1=data(savemem = 0 )
    d1.read_all_lmp( '../../8CG3/pure/trj_slab0_stack.dump' )
    d1.allframes = d1.allframes[ -200 :  : 2 ]
    df, _, stats =  violin_df_prep(d1 )
    df_8CG3 = violin_plot(df, None , sys_name='8CG3') #'seaborn_violin0.png'

    # 7KWZ Nterm
    d1=data(savemem = 0 )
    d1.read_all_lmp( '/home/xubo/Documents/test_hoomd/uconn_hpc/7KWZ_Nterm_monomer/pure_v3/trj_slab0_stack.dump' )
    d1.allframes = d1.allframes[ -200 :  : 2 ]
    df, _, stats =  violin_df_prep(d1 )
    df_Nterm = violin_plot(df, None , sys_name='7KWZ_N') #'seaborn_violin0.png'

    # Cterm
    d1=data(savemem = 0 )
    d1.read_all_lmp( '/home/xubo/Documents/test_hoomd/uconn_hpc/7KWZ_Cterm_monomer/pure_v2/trj_slab1_stack.dump' )
    d1.allframes = d1.allframes[ -200 :  : 2 ]
    df, _, stats =  violin_df_prep(d1 )
    df_Cterm = violin_plot(df, None , sys_name='7KWZ_C') #'seaborn_violin0.png'

    global_df = df_7PY2
    global_df = pd.concat([global_df, df_8CG3 ], ignore_index=True, axis= 0 )
    global_df = pd.concat([global_df, df_Nterm ], ignore_index=True, axis= 0 )
    global_df = pd.concat([global_df, df_Cterm ], ignore_index=True, axis= 0 )
    #print(global_df)

    ###################################
    ### plot global comparison    
    ###################################
    fig, ax = plt.subplots( figsize = (16, 12))
    sns.set_theme(style='whitegrid' )

    L_sys = ['7PY2', '8CG3', '7KWZ_N', '7KWZ_C' ]
    L_color = ['orange', 'pink', 'green', 'cyan' ]
    L_width = [ 0.5 , 0.5 , 0.5 , 0.5 ]
    L_bw_adjust = [ 0.5 , 0.5 , 0.5 , 0.8 ]
    for idx, sys in enumerate(L_sys):
        df = global_df[ global_df['sys'] == sys ] 
        sns.violinplot( 
            data = df,
            x = df.columns[0], 
            y = df.columns[1],
            color = L_color[idx] ,
            inner = 'box', 
            width = L_width[idx] ,
            bw_adjust = L_bw_adjust[idx],
            cut = 0 ,
            ax = ax 
        )
    
    # set axis 
    min_y = 0.8
    max_y = int ( df['N'].max() )
    ax.set_yticks(np.arange(0, max_y + 1, 5) )
    ax.set_xlabel('t', fontsize=48, labelpad=10)   
    ax.set_ylabel('N', fontsize=48, labelpad=10)   
    ax.tick_params(axis = 'both', which = 'major', labelsize = 36)


    # transparent background
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)

    plt.tight_layout()

    plt.savefig('compare.png', dpi=450, transparent= True )
    plt.close()

    stop  = timer.perf_counter()
    print('time used in sec: ', stop - start )



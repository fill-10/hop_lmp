from class_data import data
if __name__ == '__main__':
    import time as timer
    import numpy as np
    import pandas as pd 
    from class_oneframe import oneframe
    from read_1_frame import *
    ##--- trj file ---
    grofilename = '/global/scratch/xluo2/test_charmm/Ndc10-Nte10/H-Ndc-Nte-NH2/inf12/finelongtrj.gro'
    ##--- output general settings ---
    fn_prefix = '../inf12/2pt/'
    fout = fn_prefix+'nvt2pt.lammps'
    ##--- timer start ---
    start = timer.perf_counter()
    ##--- type conversion dictionary ---
    ATdict = {'N':1, 'HT1':2,'CK':5, 'C':7, 'O':8, 'NT':10, 'OW':11, 'HW1':12, 'HW2':12 }
    for key in ['CA', 'CA5', 'CB', 'CG', 'CD', 'CE', 'CZ', 'CH', 'CT', 'CI']:
        ATdict[key] = 3
    for key in ['HA1', 'HA2', 'HA3', 'HA4', 'HB1','HB2', 'HG1', 'HG2', 'HD1', 'HD2', 'HE1', 'HE2','HZ1', 'HZ2', 'HH1', 'HH2', 'HT5', 'HT6', 'HI1', 'HI2']:
        ATdict[key] = 4
    for key in ['HK1', 'HK2', 'HK3' ]:
        ATdict[key] = 6
    for key in ['OG', 'OZ', 'OI' ]:
        ATdict[key] = 6
    ##--- read in gro traj ---
    fi = open(grofilename, 'r')
    ##--- delete old ---
    import os
    try :
        os.remove(fout)
    except OSError as e:
        print(e)
    ##--- loop for all frames ---
    while 1:
        try :
            ##--- initialize a frame ---
            onef = oneframe()
            ctime, cNatom, cbox, ccols, catoms, cpos = read_1_gro(fi)
            ctime = int(round(ctime*1000/2))  # 1 ps = 1000 fs, gmx timestep = 2 fs
            onef.load_snap(ctime, cbox, catoms, ccols)
            ##--- change type to number ---
            onef.L_atom['type'].replace(ATdict, inplace=True)
            ##-- unit conversion: nm to A, ps to fs --
            onef.L_atom['xu'] = round( onef.L_atom['xu'] * 10  , 2 )
            onef.L_atom['yu'] = round( onef.L_atom['yu'] * 10  , 2 )
            onef.L_atom['zu'] = round( onef.L_atom['zu'] * 10  , 2 )
            onef.L_atom['vx'] = round( onef.L_atom['vx'] / 100 , 6 ) #1nm=10A, 1ps=1000fs
            onef.L_atom['vy'] = round( onef.L_atom['vy'] / 100 , 6 )
            onef.L_atom['vz'] = round( onef.L_atom['vz'] / 100 , 6 )
            onef.deltaX *= 10
            onef.deltaY *= 10
            onef.deltaZ *= 10
            onef.xlo *= 10
            onef.xhi *= 10
            onef.ylo *= 10
            onef.yhi *= 10
            onef.zlo *= 10
            onef.zhi *= 10
            onef.L_atom.loc[100000:, 'id'] += 100000
            ##--- save lammpstrj ---
            ##--- output file ---
            fo = open(fout, 'a')
            col_out = ['id', 'type', 'xu', 'yu', 'zu', 'vx', 'vy', 'vz']
            onef.export_lmptrj(fo, onef.L_atom, col_out)
            fo.close()
            #print(ctime, cNatom)
        except Exception as error:
            print(error)
            break
    ##--- close files ---
    fi.close()

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)


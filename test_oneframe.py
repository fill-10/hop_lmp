from class_oneframe import oneframe
from read_1_frame import *

if __name__ == '__main__':
    import time as timer
    ##--- timer starts ---
    start = timer.perf_counter()

    pdbfilename = '../nvtprod_every1ns_00-50ns.pdb'
    f = open(pdbfilename, 'r')
    for j in range(0,9):
        time, Natom, box, col, atoms, pos = read_1_pdb(f)
        onef = oneframe()
        onef.load_snap(time, box, atoms, col)
        onef.L_AN = onef.ion_gen(range(1,1+400), 1, 1, 1, [], 1, "sel[:]", 'type')
        onef.L_CT = onef.ion_gen(range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type' )
        onef.L_CT['id'] += max(onef.L_AN['id'])
        del onef.L_atom
        onef.L_atom=[]
        onef.wrap(onef.L_AN)
        onef.wrap(onef.L_CT)
        onef.find_asso(onef.L_CT, onef.L_AN, 7.8)

    time, Natom, box, col, atoms, pos = read_1_pdb(f)
    secf = oneframe()
    secf.load_snap(time, box, atoms, col)
    secf.L_AN = secf.ion_gen(range(1,1+400), 1, 1, 1, [], 1, "sel[:]", 'type')
    secf.L_CT = secf.ion_gen(range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type' )
    secf.L_CT['id'] += max(secf.L_AN['id'])
    del secf.L_atom
    secf.L_atom = []
    secf.wrap(secf.L_AN)
    secf.wrap(secf.L_CT)
    secf.find_asso(secf.L_CT, secf.L_AN,7.8, 0)
    
    f.close()
 
    ##--- timer ends ---
    end = timer.perf_counter()
    print('time used: ', end-start)
    print(secf.hist_asso_im_atom, secf.hist_asso_im_mol)
    print(secf.L_AN,secf.L_CT)
    
    """
    onef.unwrap(onef.L_AN)
    secf.unwrap(secf.L_AN)
    fastpercent = secf.findfast(secf.L_AN, onef.L_AN, 4.6)
    onef.L_AN['fast'] = secf.L_AN['fast']
    stringlength_hist = secf.findstring(secf.L_AN, onef.L_AN, 2.5)
    print('first frame: \n', onef.L_AN.loc[onef.L_AN['fast']>0,['id','mol','x','y','z','fast'] ])
    print('second frame: \n', secf.L_AN.loc[secf.L_AN['fast']>0,['id','mol', 'x','y','z','fast','string'] ])
    print(fastpercent)
    print(stringlength_hist)
    ##--- output ---
    onef.L_AN['type'] = secf.L_AN['string']
    secf.L_AN['type'] = secf.L_AN['string'] + 100

    f = open('stringsample_onef.lammpstrj','w')
    onef.export_lmptrj(f, onef.L_AN.loc[ onef.L_AN['fast']>0, ['id','mol','type','x','y','z','ix', 'iy', 'iz'] ] )
    f.close()
    f = open('stringsample_secf.lammpstrj','w')
    secf.export_lmptrj(f, secf.L_AN.loc[ secf.L_AN['fast']>0, ['id','mol','type','x','y','z','ix', 'iy', 'iz'] ] )
    f.close()
    """


def read_1_lmp(f, start=0, dim = 3):
    if start>0:
        f.seek(start)
    ###--- read controlling flags
    readtime = 0
    readnumber = 0
    readbox = 0
    readatom = 0

    item_atoms = []
    item_time = 0.0
    item_box = []
    Natom = 0
    atom_count  = 0
    col_dict = []
    # read head info
    rawline = '\n'
    while rawline !='':
        rawline = f.readline()
        cline = rawline.split()
        if cline[0] == 'ITEM:':  # control flags, find which to read
            if cline[1] == 'TIMESTEP':
                readtime = 1
                readnumber = 0
                readbox = 0
                readatom = 0
                continue
            if cline[1] == 'NUMBER':
                readnumber = 1
                readtime = 0
                readbox = 0
                readatom = 0
                continue
            if cline[1] == 'BOX':
                readbox = dim # as counter
                readtime = 0
                readnubmer = 0
                readatom = 0
                continue
            if cline[1] == 'ATOMS':
                readtime = 0
                readnumber = 0
                readbox = 0
                readatom = 1
                col_dict = cline[2:]
                if Natom>0:
                    atom_count = Natom
                    continue
                else:
                    raise Exception(f.tell(), 'Natom not given before reading atoms')
        # read info
        if readtime == 1 :
            item_time = float(cline[0])
            readtime = 0
        if readnumber == 1 :
            Natom = int(cline[0])
            atom_count = Natom
            readnumber = 0
        if readbox>0:
            # convert data type
            for k in range(0, len(cline)):
                cline[k] = float(cline[k])
            item_box += [cline]
            readbox -= 1
        if readatom and atom_count >0:
            # convert data types
            for k in range(0, len(col_dict)):
                if col_dict[k] in ['id', 'mol', 'ix', 'iy', 'iz']:
                    cline[k] = int(cline[k])
                elif col_dict[k] in ['q', 'x', 'y', 'z', 'xu', 'yu', 'zu']:
                    cline[k] = float(cline[k])
                elif col_dict[k] in ['element', 'type']:
                    try:
                        cline[k] = int(cline[k])
                    except:
                        pass
                else:
                    cline[k] = float(cline[k])
            item_atoms += [cline]
            atom_count -= 1

        if readatom and atom_count<=0:
            return item_time, Natom, item_box, col_dict, item_atoms, f.tell()

def read_1_pdb(f, start=0, dim=3):
    if start>0:
        f.seek(start)
    Natom = 0
    item_time = 0.0
    item_atoms= []
    item_box = []
    # pdb only records 'unwrapped' data
    col_dict = ['id','mol', 'type', 'xu', 'yu', 'zu']
    rawline = 'start'
    while rawline[0:6] != 'ENDMDL' and rawline !='': # EOF is ''
        rawline = f.readline()
        if rawline[0:6] == 'TITLE ':
            cline = rawline.split()
            # gromacs: time is recorded after 't='
            item_time = float( cline[cline.index('t=')+1] )
        elif rawline[0:6] == 'CRYST1':
            cline = rawline.split()
            item_box = [ [0., float(cline[1])], [0., float(cline[2])], [0, float(cline[3]) ] ]
        elif rawline[0:6] == 'ATOM  ' or rawline[0:6]=='HETATM':
            Natom += 1
            item_atoms.append(  [  int(rawline[6:11]), int(rawline[22:26]), rawline[12:16].replace(' ',''), float(rawline[30:38]), float(rawline[38:46]), float(rawline[46:54]) ] )
    
    return item_time, Natom, item_box, col_dict, item_atoms, f.tell()

def read_1_gro(f, start=0, dim=3):
    if start > 0:
        f.seek(start)
    Natom = 0
    item_time = 0.0
    item_atoms = []
    item_box   = []
    col_dict = ['mol','res', 'type','id', 'xu', 'yu', 'zu']
    rawline = 'start'
    atom_sec = 0
    n = 1
    while n>0 and rawline !='': # EOF is ''
        # if invalid .gro file, it will read to the EOF 
        # but no data collected
        rawline = f.readline()
        if atom_sec:
            item_atoms.append( \
            [   int(rawline[0:5]), \
                rawline[5:10].replace(' ',''), \
                rawline[10:15].replace(' ',''), \
                int(rawline[15:20]), \
                float(rawline[20:28]) , \
                float(rawline[28:36]) , \
                float(rawline[36:44]) , \
            ] )
            n = n-1
        else: # head
            cline = rawline.split()
            if 't=' in cline:
                item_time = float( cline[cline.index('t=')+1] )
                Natom = int( f.readline().split()[0] )
                n = Natom
                atom_sec = 1
        # box size
    cline = f.readline().split()
    item_box = [  \
                [ 0., float(cline[0]) ], \
                [ 0., float(cline[1]) ], \
                [ 0., float(cline[2]) ]  \
           ]
    #print('no box info')
    # return
    return item_time, Natom, item_box, col_dict, item_atoms, f.tell()

if __name__ == '__main__' :
    ##---test read
    """
    filename = '../test3_rest2/rest0/rest0_protein.gro'
    f = open(filename, 'r')
    while 1:
        time, Natom, box, col, atoms, pos = read_1_gro(f)
        print(time)
        print(box)
        print(col)
        print(atoms[1])
        print(atoms[-1])
        print(pos)
    f.close()
    """
    pass

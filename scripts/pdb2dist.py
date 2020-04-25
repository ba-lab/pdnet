'''
Author: Badri Adhikari, University of Missouri-St. Louis, 12-18-2019
File: Contains code to obtain distance map (pkl files) from clean PDB chains
'''

import urllib.request
from datetime import datetime
import os
import sys
import random
from math import sqrt
import numpy as np
import getopt

if sys.version_info < (3,0,0):
    print('Python 3 required!!!')
    sys.exit(1)

def usage():
    print('Usage:')
    print(sys.argv[0] + ' <-i file_pdb> <-b file_cb_dmap> <-a file_ca_dmap> [-p] [-g file_any2any_dmap(ToDo)]')

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:b:a:g:ph")
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

flag_gaps = False
file_pdb  = ''
file_cb_dmap = ''
file_ca_dmap = ''
file_aa_dmap = ''

for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-i"):
        file_pdb = os.path.abspath(a)
    elif o in ("-p"):
        flag_gaps = True
    elif o in ("-b"):
        file_cb_dmap = os.path.abspath(a)
    elif o in ("-a"):
        file_ca_dmap = os.path.abspath(a)
    elif o in ("-g"):
        file_aa_dmap = os.path.abspath(a)
    else:
        assert False, "Error!! unhandled option!!"

start = datetime.now()

if len(file_pdb) < 2:
    print('Need input pdb!')
    usage()
    sys.exit()
if len(file_cb_dmap) < 2:
    print('Need output cb dmap filename!')
    usage()
    sys.exit()

valid_amino_acids = {
    'LLP': 'K', 'TPO': 'T', 'CSS': 'C', 'OCS': 'C', 'CSO': 'C', 'PCA': 'E', 'KCX': 'K', \
    'CME': 'C', 'MLY': 'K', 'SEP': 'S', 'CSX': 'C', 'CSD': 'C', 'MSE': 'M', \
    'ALA': 'A', 'ASN': 'N', 'CYS': 'C', 'GLN': 'Q', 'HIS': 'H', 'LEU': 'L', \
    'MET': 'M', 'MHO': 'M', 'PRO': 'P', 'THR': 'T', 'TYR': 'Y', 'ARG': 'R', 'ASP': 'D', \
    'GLU': 'E', 'GLY': 'G', 'ILE': 'I', 'LYS': 'K', 'PHE': 'F', 'SER': 'S', \
    'TRP': 'W', 'VAL': 'V', 'SEC': 'U'
    }

def check_pdb_valid_row(valid_amino_acids, l):
    if (get_pdb_rname(l) in valid_amino_acids.keys()) and (l.startswith('ATOM') or l.startswith('HETA')):
        return True
    return False

def get_pdb_atom_name(l):
    return l[12: 16].strip()

def get_pdb_rnum(l):
    return int(l[22: 27].strip())

def get_pdb_rname(l):
    return l[17: 20].strip()

def get_pdb_xyz_cb(lines):
    xyz = {}
    for l in lines:
        if not (l.startswith('ATOM') or l.startswith('HE')):
            continue
        if get_pdb_atom_name(l) == 'CB':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    for l in lines:
        if not (l.startswith('ATOM') or l.startswith('HE')):
            continue
        if (get_pdb_rnum(l) not in xyz) and get_pdb_atom_name(l) == 'CA':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    return xyz

def get_pdb_xyz_ca(lines):
    xyz = {}
    for l in lines:
        if not (l.startswith('ATOM') or l.startswith('HE')):
            continue
        if get_pdb_atom_name(l) == 'CA':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    return xyz

def get_dist_maps(valid_amino_acids, file_pdb, flag_gaps, flag_any2any = False):
    f = open(file_pdb, mode = 'r')
    lines = f.read()
    f.close()
    lines = lines.splitlines()
    rnum_rnames = {}
    for l in lines:
        atom = get_pdb_atom_name(l)
        if atom != 'CA':
            continue
        if not get_pdb_rname(l) in valid_amino_acids.keys():
            print ('' + get_pdb_rname(l) + ' is unknown amino acid in ' + l)
            sys.exit(1)
        rnum_rnames[int(get_pdb_rnum(l))] = valid_amino_acids[get_pdb_rname(l)]
    seq = ''
    for i in range(max(rnum_rnames.keys())):
        if i+1 not in rnum_rnames:
            print (rnum_rnames)
            print ('Warning! ' + file_pdb + ' ! residue not defined for rnum = ' + str(i+1))
            seq = seq + '-'
        else:
            seq = seq + rnum_rnames[i+1]
    L = len(seq)
    xyz_cb = get_pdb_xyz_cb(lines)

    if not flag_gaps:
        if len(xyz_cb) != L:
            print(rnum_rnames)
            for i in range(L):
                if i+1 not in xyz_cb:
                    print('XYZ not defined for ' + str(i+1))
            print ('Error! ' + file_pdb + ' Something went wrong - len of cbxyz != seqlen!! ' + str(len(xyz_cb)) + ' ' +  str(L))
            sys.exit(1)

    cb_map = np.full((L, L), np.nan)
    for r1 in sorted(xyz_cb):
        (a, b, c) = xyz_cb[r1]
        for r2 in sorted(xyz_cb):
            (p, q, r) = xyz_cb[r2]
            cb_map[r1 - 1, r2 - 1] = sqrt((a-p)**2+(b-q)**2+(c-r)**2)
    xyz_ca = get_pdb_xyz_ca(lines)

    if not flag_gaps:
        if len(xyz_ca) != L:
            print ('Something went wrong - len of cbxyz != seqlen!! ' + str(len(xyz_ca)) + ' ' +  str(L))
            sys.exit(1)

    ca_map = np.full((L, L), np.nan)
    for r1 in sorted(xyz_ca):
        (a, b, c) = xyz_ca[r1]
        for r2 in sorted(xyz_ca):
            (p, q, r) = xyz_ca[r2]
            ca_map[r1 - 1, r2 - 1] = sqrt((a-p)**2+(b-q)**2+(c-r)**2)

    if flag_any2any:
        any_map = np.full((L, L), np.nan)
        for l1 in lines:
            if not check_pdb_valid_row(check_pdb_valid_row, l1):
                continue
            r1 = get_pdb_rnum(l1)
            (a, b, c) = (float(l1[30:38].strip()), float(l1[38:46].strip()), float(l1[46:54].strip()))
            for l2 in lines:
                if not check_pdb_valid_row(check_pdb_valid_row, l2):
                    continue
                r2 = get_pdb_rnum(l2)
                (p, q, r) = (float(l2[30:38].strip()), float(l2[38:46].strip()), float(l2[46:54].strip()))
                d = sqrt((a-p)**2+(b-q)**2+(c-r)**2)
                if any_map[r1 - 1, r2 - 1] > d:
                    any_map[r1 - 1, r2 - 1] = d
        return L, seq, cb_map, ca_map, any_map
    return L, seq, cb_map, ca_map

def parse_pdb_row(row, param):
    result = ''
    if param == 'rnum':
        result = row[22:27]
        if any(c.isalpha() for c in result):
            valid_list = {
                '   0A': '-999',
                '   0B': '-998',
                '   0C': '-997',
                '   0D': '-996',
                '   0E': '-995',
                '   0F': '-994',
                '   0G': '-993',
                '   0H': '-992',
                '   0I': '-991',
                '   0J': '-990',
            }
            if result in valid_list.keys():
                result = valid_list[result]
            else:
                #print('Alternative ' + row + ' - skipping..\n')
                return 'NA'
    if param == 'aname':
        result = row[12:16]
    if param == 'altloc':
        result = row[16:17]
    if param == 'rname':
        result = row[17:20]
    if param == 'chain':
        result = row[21:22]
    if param == 'x':
        result = row[30:38]
    if param == 'y':
        result = row[38:46]
    if param == 'z':
        result = row[46:54]
    if len(result) < 1:
        print('Error! Undefined param/result!')
        sys.exit(1)
    return result.strip()

# Input is a single PDB chain, reindexed, sorted (all clean)
print (' Extract distance maps..')

(L, seq, cb_map, ca_map) = get_dist_maps(valid_amino_acids, file_pdb, flag_gaps, False)
print(' ' + str(L))
print(' ' + str(ca_map.shape))

print (' Verify the symmetricity in distance maps..')
if not np.allclose(cb_map, cb_map.T, equal_nan = True):
    print('ERROR in symmetricity!!')
    print('Quitting..')
    sys.exit(1)

print (' Writing distance map file.. ')
np.save(file_cb_dmap, ((L, seq, cb_map.astype(np.float16))))
np.save(file_ca_dmap, ((L, seq, ca_map.astype(np.float16))))

print ('Runtime = ' + str(datetime.now() - start))

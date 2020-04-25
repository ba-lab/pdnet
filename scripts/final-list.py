'''
Author: Badri Adhikari, University of Missouri-St. Louis,  12-30-2019
File: List the IDs that match fasta, chain, and features
'''

from datetime import datetime
import os
import sys
import numpy as np
import pickle

if sys.version_info < (3,0,0):
    print('Python 3 required!!!')
    sys.exit(1)


start = datetime.now()

dir_pdb   = sys.argv[1] # './pdb/'
dir_dmap  = sys.argv[2] # './distance-map/'
dir_fasta = sys.argv[3] # './fasta/'
dir_pkl   = sys.argv[4] # './features/'

valid_amino_acids = {
    'LLP': 'K', 'TPO': 'T', 'CSS': 'C', 'OCS': 'C', 'CSO': 'C', 'PCA': 'E', 'KCX': 'K', \
    'CME': 'C', 'MLY': 'K', 'SEP': 'S', 'CSX': 'C', 'CSD': 'C', 'MSE': 'M', \
    'ALA': 'A', 'ASN': 'N', 'CYS': 'C', 'GLN': 'Q', 'HIS': 'H', 'LEU': 'L', \
    'MET': 'M', 'MHO': 'M', 'PRO': 'P', 'THR': 'T', 'TYR': 'Y', 'ARG': 'R', 'ASP': 'D', \
    'GLU': 'E', 'GLY': 'G', 'ILE': 'I', 'LYS': 'K', 'PHE': 'F', 'SER': 'S', \
    'TRP': 'W', 'VAL': 'V', 'SEC': 'U'
    }

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

def chain2sequence(valid_amino_acids, chain):
    f = open(chain, mode = 'r')
    lines = f.read()
    f.close()
    lines = lines.splitlines()

    end = -1
    for line in lines:
        if len(line) < 10:
            continue
        if parse_pdb_row(line, 'aname') != 'CA':
            continue
        if parse_pdb_row(line, 'rname') not in valid_amino_acids.keys():
            print('ERROR! ' + parse_pdb_row(line, 'rname') + ' residue not defined!')
            sys.exit(1)
        end = int(parse_pdb_row(line, 'rnum'))

    seq = ''
    for line in lines:
        if len(line) < 10:
            continue
        if parse_pdb_row(line, 'aname') != 'CA':
            continue
        if parse_pdb_row(line, 'rname') not in valid_amino_acids.keys():
            print('ERROR! ' + parse_pdb_row(line, 'rname') + ' residue not defined!')
            sys.exit(1)
        seq += valid_amino_acids[parse_pdb_row(line ,'rname')]

    if len(seq) != end:
        print('ERROR! len(seq) ' + str(len(seq)) + ' is unequal to max residue number ' + str(end))
        sys.exit(1)

    if len(seq) < 3:
        print('WARNING!! Too short sequence')

    return seq

for filename in sorted(os.listdir(dir_fasta)):
    if not filename.endswith(".fasta"):
        continue
    id = filename.split('.')[0]

    seq1 = chain2sequence(valid_amino_acids, dir_pdb + id + '.pdb')
    seq2 = ''
    with open(dir_fasta + id + '.fasta') as ffasta:
        for fastaline in ffasta:
            if (fastaline.startswith('>')):
                continue
            seq2 = fastaline.strip()
    if seq1 != seq2:
        print(id)
        print(len(seq1), len(seq2))
        print(seq1)
        print(seq2)
        print('ERROR! chain seq mismatches fasta seq!!')
        sys.exit(1)

    (L, seq, cb_map) = np.load(dir_dmap + id + '-cb.npy', allow_pickle = True)
    if seq1 != seq:
        print(id)
        print(len(seq1), len(seq))
        print(seq1)
        print(seq)
        print('ERROR! chain seq mismatches dmap seq!!')
        sys.exit(1)

    if len(sys.argv) > 4:
        features = pickle.load(open(dir_pkl + id + '.pkl', 'rb'))
        l = len(features['seq'])
        seq3 = features['seq']
        if seq1 != seq3:
            print(id)
            print(seq1)
            print(seq3)
            print(len(seq1), len(seq3))
            print('ERROR! chain seq mismatches feature sequence!!')
            sys.exit(1)

    print(id, L)

'''
Author: Badri Adhikari, University of Missouri-St. Louis, 12-30-2019
File: Contains code to download and clean PDBs and prepare fasta and distance maps
Input: PDB ID list and working directories
Outputs: 1) Reindexed clean PDB files
         2) Corresponding fasta files
         3) Numpy distance maps
'''

import urllib.request
from datetime import datetime
import os
import sys
import random
from math import sqrt
import numpy as np

if sys.version_info < (3,0,0):
    print('Python 3 required!!!')
    sys.exit(1)

if len(sys.argv) < 8:
    print('Error! Insufficient parameters!')
    sys.exit(1)

start = datetime.now()

input_list         = sys.argv[1] # './pdnet.lst'
dir_pdb            = sys.argv[2] # './pdb/'
dir_dmap           = sys.argv[3] # './distance-map/'
dir_fasta          = sys.argv[4] # './fasta/'
dir_downloads      = sys.argv[5] # './temp-downloads/'
# pisces db has 1000+ sequences with 512+ residues
max_chain_size     = int(sys.argv[6]) # 512
# pisces db has all 50+ but cath as as low as 14
min_chain_size     = int(sys.argv[7]) # 12

valid_amino_acids = {
    'LLP': 'K', 'TPO': 'T', 'CSS': 'C', 'OCS': 'C', 'CSO': 'C', 'PCA': 'E', 'KCX': 'K', \
    'CME': 'C', 'MLY': 'K', 'SEP': 'S', 'CSX': 'C', 'CSD': 'C', 'MSE': 'M', \
    'ALA': 'A', 'ASN': 'N', 'CYS': 'C', 'GLN': 'Q', 'HIS': 'H', 'LEU': 'L', \
    'MET': 'M', 'MHO': 'M', 'PRO': 'P', 'THR': 'T', 'TYR': 'Y', 'ARG': 'R', 'ASP': 'D', \
    'GLU': 'E', 'GLY': 'G', 'ILE': 'I', 'LYS': 'K', 'PHE': 'F', 'SER': 'S', \
    'TRP': 'W', 'VAL': 'V', 'SEC': 'U'
    }
    #, 'UNK': 'X', 'PYL': 'X'

flag_check_for_structural_gaps = False

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

def get_dist_maps(valid_amino_acids, file_pdb):
    f = open(file_pdb, mode = 'r')
    lines = f.read()
    f.close()
    lines = lines.splitlines()
    rnum_rnames = {}
    for l in lines:
        atom = get_pdb_atom_name(l)
        if atom != 'CA':
            continue
        #if int(get_pdb_rnum(l)) in rnum_rnames:
            #warnings.warn ('Warning!! ' + file_pdb + ' - multiple CA rows - rnum = ' + str(get_pdb_rnum(l)))
        if not get_pdb_rname(l) in valid_amino_acids.keys():
            print ('' + get_pdb_rname(l) + ' is unknown amino acid in ' + l)
            sys.exit(1)
        rnum_rnames[int(get_pdb_rnum(l))] = valid_amino_acids[get_pdb_rname(l)]
    seq = ""
    for i in range(max(rnum_rnames.keys())):
        if i+1 not in rnum_rnames:
            print (rnum_rnames)
            print ('Error! ' + file_pdb + ' ! residue not defined for rnum = ' + str(i+1))
            sys.exit (1)
        seq = seq + rnum_rnames[i+1]
    L = len(seq)
    xyz_cb = get_pdb_xyz_cb(lines)
    if len(xyz_cb) != L:
        print(rnum_rnames)
        for i in range(L):
            if i+1 not in xyz_cb:
                print('XYZ not defined for ' + str(i+1))
        print ('Error! ' + file_pdb + ' Something went wrong - len of cbxyz != seqlen!! ' + str(len(xyz_cb)) + ' ' +  str(L))
        sys.exit(1)
    cb_map = np.zeros((L, L))
    for r1 in sorted(xyz_cb):
        (a, b, c) = xyz_cb[r1]
        for r2 in sorted(xyz_cb):
            (p, q, r) = xyz_cb[r2]
            cb_map[r1 - 1, r2 - 1] = sqrt((a-p)**2+(b-q)**2+(c-r)**2)
    xyz_ca = get_pdb_xyz_ca(lines)
    if len(xyz_ca) != L:
        print ('Something went wrong - len of cbxyz != seqlen!! ' + str(len(xyz_ca)) + ' ' +  str(L))
        sys.exit(1)
    ca_map = np.zeros((L, L))
    for r1 in sorted(xyz_ca):
        (a, b, c) = xyz_ca[r1]
        for r2 in sorted(xyz_ca):
            (p, q, r) = xyz_ca[r2]
            ca_map[r1 - 1, r2 - 1] = sqrt((a-p)**2+(b-q)**2+(c-r)**2)
    '''
    any_map = np.full((L, L), np.inf)
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
    '''
    return L, seq, cb_map, ca_map

# Reindex Chain. alternative locations removed, one model, one chain, sudden decrease in residue number trimmed
def reindex_chain_from_rnum1(valid_amino_acids, inpdb, maxresnum = 512):
    res_fiff = 0;
    atom_counter = 0;
    prev_res_num_in_inPDB = -10000;
    this_res_num_in_inPDB = 0;

    f = open(inpdb, mode = 'r')
    lines = f.read()
    f.close()
    lines = lines.splitlines()
    print('  input rows:', len(lines))

    # Skip the last lines of HETATM
    skip_from_this_line_on = '--';
    for line in reversed(lines):
        if line.startswith('HETATM') and 'MSE' not in line:
            skip_from_this_line_on = line
        else:
            break

    print('  skipping lines after: ', skip_from_this_line_on)
    new_rnum = 0;
    # skip all the residues that do not have CA
    residues_to_skip = {}
    for line in lines:
        if len(line) < 10:
            continue
        if not (line[16:17] == ' ' or line[16:17] == 'A'):
            continue
        this_rnum = parse_pdb_row(line, 'rnum')
        residues_to_skip[this_rnum] = 1
    for line in lines:
        if len(line) < 10:
            continue
        if not (line[16:17] == ' ' or line[16:17] == 'A'):
            continue
        this_rnum = parse_pdb_row(line, 'rnum')
        this_aname = parse_pdb_row(line, 'aname')
        if this_aname == 'CA' and this_rnum in residues_to_skip:
            del residues_to_skip[this_rnum]

    lines_to_write = []
    for line in lines:
        if len(line) < 10:
            continue
        if not (line[16:17] == ' ' or line[16:17] == 'A'):
            continue
        if line == skip_from_this_line_on:
            break
        this_rnum = parse_pdb_row(line, 'rnum')
        if this_rnum in residues_to_skip:
            continue
        this_rname = parse_pdb_row(line, 'rname')
        if this_rname not in valid_amino_acids.keys():
            continue
        if this_rnum == 'NA':
            continue
        if prev_res_num_in_inPDB != this_rnum:
            prev_res_num_in_inPDB = this_rnum
            new_rnum += 1
        if new_rnum == maxresnum + 1:
            break
        atom_counter += 1
        str_atom_counter = '%5s' % (atom_counter)
        str_new_rnum = '%4s' % (new_rnum)
        lines_to_write.append(line[:6] + str_atom_counter + line[11:16] + ' ' + line[17:20] + '  ' + str_new_rnum + ' ' + line[27:])

    print('  output rows:', len(lines_to_write))

    if (len(lines_to_write) < 10):
        print('WARNING! Too few lines to write.. \n')

    f = open(inpdb + '.tmp', mode = 'w')
    for line in lines_to_write:
        f.write(line + '\n')
    f.close()

    os.system('mv' + ' ' + inpdb + '.tmp' + ' ' + inpdb)
    return new_rnum

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


f = open(input_list, mode = 'r')
input_pdbs = f.read()
f.close()
input_pdbs = input_pdbs.splitlines()

for ipdb in range(len(input_pdbs)):
    sys.stdout.flush()
    if input_pdbs[ipdb].startswith('IDs'):
        continue
    #if ipdb == 5:
    #    continue
    cols = input_pdbs[ipdb].strip()
    id = cols[0:4]
    chain = cols[4]

    '''
    # Cross checking
    if os.path.isfile(dir_fasta + id + chain + '.fasta'):
        continue
    '''

    print('')
    print(ipdb, '/', len(input_pdbs))
    print(id, chain)
    #'''
    # Skip if already done
    if os.path.isfile(dir_fasta + id + chain + '.fasta') and os.path.isfile(dir_pdb + id + chain + '.pdb'):
        seq1 = chain2sequence(valid_amino_acids, dir_pdb + id + chain + '.pdb')
        if len(seq1) > min_chain_size:
            seq2 = ''
            with open(dir_fasta + id + chain + '.fasta') as ffasta:
                for fastaline in ffasta:
                    if (fastaline.startswith('>')):
                        continue
                    seq2 = fastaline.strip()
            if seq1 == seq2:
                print ('Already done.. skipping..')
                continue

    if (not os.path.isfile(dir_downloads + id + '.pdb')):
        print (' Downloading..')
        try:
            urllib.request.urlretrieve("https://files.rcsb.org/download/" + id + '.pdb', dir_downloads + id + '.pdb')
        except:
            print('Something went wrong! Could not download ' + id + ' at ' + "https://files.rcsb.org/download/" + id + '.pdb')
            continue

    print (' Extract chain..')
    fchain = open(dir_pdb + id + chain + '.pdb', "w")
    chain_char = chain
    if chain_char == '0':
        chain_char = ' '

    with open(dir_downloads + id + '.pdb') as fpdb:
        for atomline in fpdb:
            if (not (atomline.startswith('ATOM') or atomline.startswith('HETATM'))):
                continue
            if (len(atomline) < 22):
                continue
            if (atomline[21] != chain_char):
                continue
            fchain.write(atomline)
    fchain.close()

    print (' Reindex chain..')
    total_residues = reindex_chain_from_rnum1(valid_amino_acids, dir_pdb + id + chain + '.pdb', max_chain_size)
    if total_residues < min_chain_size:
        print('WARNING! Too few residues.. skipping')
        continue
    #'''
    id = id + chain

    print (' Extract distance maps..')
    (L, seq, cb_map, ca_map) = get_dist_maps(valid_amino_acids, dir_pdb + id + '.pdb')
    print(' ' + str(L))
    print(' ' + str(ca_map.shape))

    if flag_check_for_structural_gaps:
        print (' Check for large structural gaps..')
        for i in range(L-1):
            if (ca_map[i, i+1] > 6.0 or ca_map[i, i+1] < 1.0 ):
                print (' WARNING! distance between residue index ' + str(i+1) +  ' and ' + str(i+2) + ' is ' + str(ca_map[i, i+1]) + ' Angstroms')
                print (' Trimming pdb..')
                print (' (again) Reindex chain..')
                reindex_chain_from_rnum1(valid_amino_acids, dir_pdb + id + '.pdb', i+1)
                print (' (again) Extract distance maps..')
                (L, seq, cb_map, ca_map) = get_dist_maps(valid_amino_acids, dir_pdb + id + '.pdb')
                print(' ' + str(L))
                break

    print (' Write fasta file..')
    f = open(dir_fasta + id + '.fasta', mode = 'w')
    f.write('>' + id + '\n')
    f.write(chain2sequence(valid_amino_acids, dir_pdb + id + '.pdb') + '\n')
    f.close()

    print (' Verify the fasta sequence..')
    if (not os.path.isfile(dir_fasta + id + '.fasta')):
        print(' ERROR: ' + id + '.fasta not found!')
    sequence = ''
    with open(dir_fasta + id + '.fasta') as ffasta:
        for fastaline in ffasta:
            if (fastaline.startswith('>')):
                continue
            sequence = fastaline.strip()

    print (' Verify distance map length and sequence length .. ')
    if L != len(sequence):
        print('ERROR', 'L is', L, 'seq-len is', len(sequence))

    print (' Verify the symmetricity in distance maps..')
    for p in range(0, L):
        for q in range(0, L):
            if cb_map[p, q] != cb_map[q, p]:
                print('ERROR in symmetricity!!')
                sys.exit(1)

    print (' Writing distance map file.. ')
    np.save(dir_dmap + id + '-cb', ((L, seq, cb_map.astype(np.float16))))
    np.save(dir_dmap + id + '-ca', ((L, seq, ca_map.astype(np.float16))))

print ('Runtime = ' + str(datetime.now() - start))

'''
Author: Badri Adhikari, University of Missouri-St. Louis,  12-18-2019
File: Generate features for fasta that do not match the PDB fasta
      i.e. reuse existing features for most proteins
'''

from datetime import datetime
import os
import sys
import numpy as np
import pickle

if sys.version_info < (3,0,0):
    print('Python 3 required!!!')
    sys.exit(1)

input_list            = sys.argv[1]
dir_fasta             = sys.argv[2]
dir_precomputed_feats = sys.argv[3]
dir_precomputed_pssm  = sys.argv[4]
dir_out_jobs          = sys.argv[5]
dir_out_feats         = sys.argv[6]

f = open(input_list, 'r')
for line in f.readlines():
    id = line.strip().split()[0]
    print(id)
    os.system('mkdir -p ' + dir_out_jobs + id)
    os.system('touch ' + dir_out_jobs + id + '/' + id + '.a3m' )
    os.system('cp ' + dir_precomputed_feats + id + '/' + 'psipred/' + id + '.ss2' +                ' ' + dir_out_jobs + id + '/' + id + '.ss2' )
    os.system('cp ' + dir_precomputed_feats + id + '/' + 'psipred/' + id + '.solv' +               ' ' + dir_out_jobs + id + '/' + id + '.solv' )
    os.system('cp ' + dir_precomputed_feats + id + '/' + 'alnstat/' + id + '.aln' +                ' ' + dir_out_jobs + id + '/' + id + '.aln' )
    os.system('cp ' + dir_precomputed_feats + id + '/' + 'alnstat/' + id + '.colstats' +           ' ' + dir_out_jobs + id + '/' + id + '.colstats' )
    os.system('cp ' + dir_precomputed_feats + id + '/' + 'alnstat/' + id + '.pairstats' +          ' ' + dir_out_jobs + id + '/' + id + '.pairstats' )
    os.system('cp ' + dir_precomputed_feats + id + '/' + 'ccmpred/' + id + '.ccmpred' +            ' ' + dir_out_jobs + id + '/' + id + '.ccmpred' )
    os.system('cp ' + dir_precomputed_feats + id + '/' + 'freecontact/' + id + '.freecontact.rr' + ' ' + dir_out_jobs + id + '/' + id + '.freecontact.rr' )
    sys.stdout.flush()
    os.system('python3 ./scripts/make-features.py -f ' + dir_fasta + id + '.fasta' + ' -j ' + dir_out_jobs + id + ' -p ' + dir_out_feats + id + '.pkl' + ' ' + ' > ' + dir_out_jobs + id + '.log')

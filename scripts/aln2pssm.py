#!/usr/bin/python
# Badri Adhikari, 3/1/2020

import sys
import numpy as np
import string
import tensorflow as tf
import os
import random

''' See below as well
'''
tf.compat.v1.enable_eager_execution()
#tf.set_random_seed(1)

os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

inaln = sys.argv[1]
outpssm = sys.argv[2]

# read A3M and convert letters into
# integers in the 0..20 range
def parse_a3m(filename):
    seqs = []
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    # read file line by line
    for line in open(filename,"r"):
        # skip labels
        if line[0] != '>':
            # remove lowercase letters and right whitespaces
            seqs.append(line.rstrip().translate(table))
    # convert letters into numbers
    alphabet = np.array(list("ARNDCQEGHILKMFPSTWYV-"), dtype='|S1').view(np.uint8)
    msa = np.array([list(s) for s in seqs], dtype='|S1').view(np.uint8)
    for i in range(alphabet.shape[0]):
        msa[msa == alphabet[i]] = i
    # treat all unknown characters as gaps
    msa[msa > 20] = 20
    return msa

# 1-hot MSA to PSSM
def msa2pssm(msa1hot, w):
    beff = tf.reduce_sum(w)
    f_i = tf.reduce_sum(w[:,None,None]*msa1hot, axis=0) / beff + 1e-9
    h_i = tf.reduce_sum( -f_i * tf.math.log(f_i), axis=1)
    return tf.concat([f_i, h_i[:,None]], axis=1)

# reweight MSA based on cutoff
def reweight(msa1hot, cutoff):
    with tf.name_scope('reweight'):
        id_min = tf.cast(tf.shape(msa1hot)[1], tf.float32) * cutoff
        id_mtx = tf.tensordot(msa1hot, msa1hot, [[1,2], [1,2]])
        id_mask = id_mtx > id_min
        w = 1.0/tf.reduce_sum(tf.cast(id_mask, dtype=tf.float32),-1)
    return w

a3m = parse_a3m(inaln)
msa1hot  = tf.one_hot(a3m, 21, dtype=tf.float32)
w = reweight(msa1hot, 0.8)
print('msa1hot', msa1hot.shape, 'w', w.shape)

f1d_pssm = msa2pssm(msa1hot, w)
print('f1d_pssm', f1d_pssm.shape)

ncol = a3m.shape[1]
nrow = a3m.shape[0]

f1d = tf.concat(values=[f1d_pssm], axis=1)
f1d = tf.expand_dims(f1d, axis=0)
f1d = tf.reshape(f1d, [ncol,22])

''' May need to comment this for enable_eager_execution
sess = tf.Session()
with sess.as_default():
    f1d = f1d.eval()
'''

np.save(outpssm, f1d.numpy().astype(np.float16))

print(f"saved {outpssm}")

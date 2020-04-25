'''
Author: Badri Adhikari, University of Missouri-St. Louis,  1-25-2020
File: Contains the code to train and test learning real-valued distances, binned-distances and contact maps
'''

import tensorflow as tf
from tensorflow.keras.callbacks import ModelCheckpoint
import os
import sys
import numpy as np
import datetime
import argparse
#%matplotlib inline

sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])) + '/lib')

from dataio import *
from metrics import *
from generator import *
from models import *
from losses import *

flag_plots = False

if flag_plots:
    from plots import *

if sys.version_info < (3,0,0):
    print('Python 3 required!!!')
    sys.exit(1)

# Allow GPU memory growth
if hasattr(tf, 'GPUOptions'):
    import keras.backend as K
    gpu_options = tf.GPUOptions(allow_growth=True)
    sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))
    K.tensorflow_backend.set_session(sess)
else:
    # For other GPUs
    for gpu in tf.config.experimental.list_physical_devices('GPU'):
        tf.config.experimental.set_memory_growth(gpu, True)

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', type=str, required = True, dest = 'type', help="contact/distance/binned")
    parser.add_argument('-w', type=str, required = True, dest = 'file_weights', help="hdf5 weights file")
    parser.add_argument('-n', type=int, required = True, dest = 'num_chains_for_training', help="number of pdbs to use for training (use -1 for ALL)")
    parser.add_argument('-c', type=int, required = True, dest = 'training_window', help="crop size (window) for training, 64, 128, etc. ")
    parser.add_argument('-e', type=int, required = True, dest = 'training_epochs', help="# of epochs")
    parser.add_argument('-d', type=int, required = True, dest = 'arch_depth', help="residual arch depth")
    parser.add_argument('-f', type=int, required = True, dest = 'filters_per_layer', help="number of convolutional filters in each layer")
    parser.add_argument('-p', type=str, required = True, dest = 'dir_dataset', help="path where all the data (including .lst) is located")
    parser.add_argument('-v', type=int, required = True, dest = 'flag_evaluate_only', help="1 = don't train only evaluate")
    parser.add_argument('-o', type=int, required = True, dest = 'flag_evaluate_cameo', help="1 = evaluate on CAMEO-hard dataset as well (takes long)")
    args = parser.parse_args()
    return args

args = get_args()

contact_or_dist_or_binned = args.type
file_weights              = args.file_weights
num_chains_for_training   = args.num_chains_for_training
training_window           = args.training_window
training_epochs           = args.training_epochs
arch_depth                = args.arch_depth
filters_per_layer         = args.filters_per_layer
dir_dataset               = args.dir_dataset
flag_evaluate_only        = False
if args.flag_evaluate_only == 1:
    flag_evaluate_only = True
flag_evaluate_cameo       = False
if args.flag_evaluate_cameo == 1:
    flag_evaluate_cameo = True
pad_size                  = 10
batch_size                = 2
expected_n_channels       = 57

print('Start ' + str(datetime.datetime.now()))

print('')
print('Parameters:')
print('contact_or_dist_or_binned', contact_or_dist_or_binned)
print('num_chains_for_training', num_chains_for_training)
print('file_weights', file_weights)
print('training_window', training_window)
print('training_epochs', training_epochs)
print('arch_depth', arch_depth)
print('filters_per_layer', filters_per_layer)
print('pad_size', pad_size)
print('batch_size', batch_size)
print('dir_dataset', dir_dataset)

if not (contact_or_dist_or_binned == 'distance' or contact_or_dist_or_binned == 'contact' or contact_or_dist_or_binned == 'binned'):
    print('ERROR! Invalid input choice!')
    sys.exit(1)

all_feat_paths = [dir_dataset + '/deepcov/features/', dir_dataset + '/psicov/features/', dir_dataset + '/cameo/features/']
all_dist_paths = [dir_dataset + '/deepcov/distance/', dir_dataset + '/psicov/distance/', dir_dataset + '/cameo/distance/']


bins = {}
if contact_or_dist_or_binned == 'binned':
    bins[0] = '0.0 4.0'
    b = 1
    range_min = float(bins[0].split()[1])
    interval = 0.2
    while(range_min < 8.0):
        bins[b] = str(round(range_min, 2)) + ' ' + str(round(range_min + interval, 2))
        range_min += interval
        b += 1
    while(range_min <= 26):
        interval += 0.2
        bins[b] = str(round(range_min, 2)) + ' ' + str(round(range_min + interval, 2))
        b += 1
        range_min += interval
    bins[b] = '26.0 1000.0'
    print('')
    print('Number of bins', len(bins))
    print('Actual bins:', bins)

deepcov_list = load_list(dir_dataset + '/deepcov.lst', num_chains_for_training)

length_dict = {}
for pdb in deepcov_list:
    (ly, seqy, cb_map) = np.load(dir_dataset + '/deepcov/distance/' + pdb + '-cb.npy', allow_pickle = True)
    length_dict[pdb] = ly

print('')
print('Split into training and validation set..')
valid_pdbs = deepcov_list[:100]
train_pdbs = deepcov_list[100:]

print('Total validation proteins : ', len(valid_pdbs))
print('Total training proteins   : ', len(train_pdbs))

print('')
print('Validation proteins: ', valid_pdbs)

train_generator = ''
valid_generator = ''
if contact_or_dist_or_binned == 'distance':
    train_generator = DistGenerator(train_pdbs, all_feat_paths, all_dist_paths, training_window, pad_size, batch_size, expected_n_channels)
    valid_generator = DistGenerator(valid_pdbs, all_feat_paths, all_dist_paths, training_window, pad_size, batch_size, expected_n_channels)
if contact_or_dist_or_binned == 'contact':
    train_generator = ContactGenerator(train_pdbs, all_feat_paths, all_dist_paths, training_window, pad_size, batch_size, expected_n_channels)
    valid_generator = ContactGenerator(valid_pdbs, all_feat_paths, all_dist_paths, training_window, pad_size, batch_size, expected_n_channels)
if contact_or_dist_or_binned == 'binned':
    train_generator = BinnedDistGenerator(train_pdbs, all_feat_paths, all_dist_paths, bins, training_window, pad_size, batch_size, expected_n_channels)
    valid_generator = BinnedDistGenerator(valid_pdbs, all_feat_paths, all_dist_paths, bins, training_window, pad_size, batch_size, expected_n_channels)

print('')
print('len(train_generator) : ' + str(len(train_generator)))
print('len(valid_generator) : ' + str(len(valid_generator)))

X, Y = train_generator[1]
print('Actual shape of X    : ' + str(X.shape))
print('Actual shape of Y    : ' + str(Y.shape))

print('')
print('Channel summaries:')
summarize_channels(X[0, :, :, :], Y[0])

if flag_plots:
    print('')
    print('Inputs/Output of protein', 0)
    plot_protein_io(X[0, :, :, :], Y[0, :, :, 0])

print('')
print('Build a model..')
model = ''
if contact_or_dist_or_binned == 'distance':
    model = deepcon_rdd_distances(training_window, arch_depth, filters_per_layer, expected_n_channels)
if contact_or_dist_or_binned == 'contact':
    model = deepcon_rdd(training_window, arch_depth, filters_per_layer, expected_n_channels)
if contact_or_dist_or_binned == 'binned':
    model = deepcon_rdd_binned(training_window, arch_depth, filters_per_layer, len(bins), expected_n_channels)

print('')
print('Compile model..')
if contact_or_dist_or_binned == 'distance':
    model.compile(loss = 'logcosh', optimizer = 'rmsprop', metrics = ['mae'])
if contact_or_dist_or_binned == 'contact':
    model.compile(loss = 'binary_crossentropy', optimizer = 'rmsprop', metrics = ['accuracy'])
if contact_or_dist_or_binned == 'binned':
    model.compile(loss = 'categorical_crossentropy', optimizer = 'rmsprop', metrics = ['accuracy'])

print(model.summary())

if flag_evaluate_only == 0:
    if contact_or_dist_or_binned == 'binned':
        print('')
        print('Load weights from distance model..')
        distmodel = deepcon_rdd_distances(training_window, arch_depth, filters_per_layer, expected_n_channels)
        distmodel.load_weights('distance256.hdf5')
        for i, layer in enumerate(distmodel.layers):
            print(' Copying weights of layer ', i, layer.name)
            model.layers[i].set_weights(layer.get_weights())
            if layer.name == 'activation_516':
                break
    if os.path.exists(file_weights):
        print('')
        print('Loading existing weights..')
        model.load_weights(file_weights)

    print('')
    print('Train..')

    history = model.fit_generator(generator = train_generator,
        validation_data = valid_generator,
        callbacks = [ModelCheckpoint(filepath = file_weights, monitor = 'val_loss', save_best_only = True, save_weights_only = True, verbose = 1)],
        verbose = 1,
        max_queue_size = 8,
        workers = 1,
        use_multiprocessing = False,
        shuffle = True ,
        epochs = training_epochs)

    if flag_plots:
        plot_learning_curves(history)

LMAX = 512

print('')
print('Evaluate validation set..')

if contact_or_dist_or_binned == 'distance':
    model = deepcon_rdd_distances(LMAX, arch_depth, filters_per_layer, expected_n_channels)
    model.load_weights(file_weights)
    eval_distance_predictions(model, valid_pdbs, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, False, LMAX, expected_n_channels)

if contact_or_dist_or_binned == 'contact':
    model = deepcon_rdd(LMAX, arch_depth, filters_per_layer, expected_n_channels)
    model.load_weights(file_weights)
    eval_contact_predictions(model, valid_pdbs, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, False, LMAX, expected_n_channels)

if contact_or_dist_or_binned == 'binned':
    model = deepcon_rdd_binned(LMAX, arch_depth, filters_per_layer, len(bins), expected_n_channels)
    model.load_weights(file_weights)
    eval_binned_predictions(model, valid_pdbs, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, False, LMAX, bins, expected_n_channels)

print('')
print('Evaluate test set..')

psicov_list = load_list(dir_dataset + 'psicov.lst')

length_dict = {}
for pdb in psicov_list:
    (ly, seqy, cb_map) = np.load(dir_dataset + '/psicov/distance/' + pdb + '-cb.npy', allow_pickle = True)
    length_dict[pdb] = ly

if contact_or_dist_or_binned == 'distance':
    model = deepcon_rdd_distances(LMAX, arch_depth, filters_per_layer, expected_n_channels)
    model.load_weights(file_weights)
    eval_distance_predictions(model, psicov_list, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, True, LMAX, expected_n_channels)

if contact_or_dist_or_binned == 'contact':
    model = deepcon_rdd(LMAX, arch_depth, filters_per_layer, expected_n_channels)
    model.load_weights(file_weights)
    eval_contact_predictions(model, psicov_list, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, True, LMAX, expected_n_channels)

if contact_or_dist_or_binned == 'binned':
    model = deepcon_rdd_binned(LMAX, arch_depth, filters_per_layer, len(bins), expected_n_channels)
    model.load_weights(file_weights)
    eval_binned_predictions(model, psicov_list, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, True, LMAX, bins, expected_n_channels)

if flag_evaluate_cameo:
    LMAX = 1300

    print('')
    print('Evaluate on CAMEO-HARD set..')

    cameo_list = load_list(dir_dataset + 'cameo-hard.lst')

    length_dict = {}
    for pdb in cameo_list:
        (ly, seqy, cb_map) = np.load(dir_dataset + '/cameo/distance/' + pdb + '-cb.npy', allow_pickle = True)
        length_dict[pdb] = ly

    if contact_or_dist_or_binned == 'distance':
        model = deepcon_rdd_distances(LMAX, arch_depth, filters_per_layer, expected_n_channels)
        model.load_weights(file_weights)
        eval_distance_predictions(model, cameo_list, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, True, LMAX, expected_n_channels)

    if contact_or_dist_or_binned == 'contact':
        model = deepcon_rdd(LMAX, arch_depth, filters_per_layer, expected_n_channels)
        model.load_weights(file_weights)
        eval_contact_predictions(model, cameo_list, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, True, LMAX, expected_n_channels)

    if contact_or_dist_or_binned == 'binned':
        model = deepcon_rdd_binned(LMAX, arch_depth, filters_per_layer, len(bins), expected_n_channels)
        model.load_weights(file_weights)
        eval_binned_predictions(model, cameo_list, length_dict, all_feat_paths, all_dist_paths, pad_size, flag_plots, True, LMAX, bins, expected_n_channels)

print('')
print ('Everything done! ' + str(datetime.datetime.now()) )

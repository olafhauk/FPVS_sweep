"""
===========
Config file for FPVS with Frequency Sweep
===========
"""

import os
import sys

import numpy as np

###############################################################################
# paths to data.

# path to acquired raw data
cbu_path = '/megdata/cbu/fpvs'

# path to data for pre-processing
data_path = '/group/erp/data/olaf.hauk/MEG/FPVS/data_sweep'

if not os.path.isdir(data_path): # create if necessary
    os.mkdir(data_path)

# Compute TRF (1) or not (0)
do_tfr = 0

###############################################################################
# Mapping betwen filenames and subjects

map_subjects = {
    1: ('meg19_0380', '191008'),  # pilot frequency sweep
    2: ('meg19_0381', '191015'),  # pilot frequency sweep
    3: ('meg19_0383', '191018'),  # pilot frequency sweep
    4: ('meg19_0384', '191021'),  # pilot frequency sweep
}

# create subject-specific data directories if necessary
for ss in map_subjects:
    subjects_dir = os.path.join(data_path, map_subjects[ss][0]) # subject data dir
    if not os.path.isdir(subjects_dir):
        print('Creating directory %s.' % subjects_dir)
        os.mkdir(subjects_dir)
    fig_dir = os.path.join(data_path, map_subjects[ss][0], 'Figures') # subject figure dir
    if not os.path.isdir(fig_dir):
        print('Creating directory %s.' % fig_dir)
        os.mkdir(fig_dir)


###############################################################################
# Bad channels

bad_channels = {
    1: {'eeg': [],
        'meg': ['MEG1123', 'MEG2223', 'MEG0813']},
    2: {'eeg': [],
        'meg': ['MEG1123', 'MEG0093', 'MEG0813', 'MEG1412']},
    3: {'eeg': [],
        'meg': ['MEG2312', 'MEG2311', 'MEG0413']},
    4: {'eeg': ['EEG012'],
        'meg': ['MEG0723', 'MEG1123', 'MEG0813']},
}

# For subjects without clean ECG channel,
# use the following magnetometers in ICA (else specify '' to use ECG)
ECG_channels = {
    1: '',
    2: '',
    3: '',
    4: '',
}


###############################################################################
# FPVS

# time segment to remove at beginning of run (s)
fpvs_leadin = 1.

# The actual sweep frequencies (consistent with fpvs_n_sweeps)
# Note: faces only have one frequency, 6Hz, no sweeps
fpvs_freqs = [12., 10., 6., 4., 3.]
fpvs_n_sweeps = 5

# oddball frequencies
# Note: at some point in scripts frequencies will be rounded to 2 decimal
# points
fpvs_odd_freq = {}
fpvs_odd_freq['faces'] = 1.2
fpvs_odd_freq['words'] = 1.

# number of harmonics to sum up
fpvs_n_harms = 20

# duration of frequency segment per run
# fpvs_n_sweeps*fpvs_sweep_duration is the run duration
fpvs_sweep_duration = 12.

# event markers for run onsets
fpvs_event_ids = [14, 15, 16, 17]

# Frequencies for PSD plots
# # frequencies as times in seconds
# times = np.array([1.0, 1.25, 1.66, 2., 6., 7.5, 10.])/1000.
# for faces 1.2 Hz will be used instead of 1Hz
topo_times = {}
topo_times['words'] = np.array([1., 2., 3., 4., 6., 10., 12.]) / 1000.
topo_times['faces'] = np.array([1.2, 2.4, 3.6, 4.8, 6., 12.]) / 1000.

# for plot_joint() display
crop_times = [0.0005, 0.013]

# extract effects at target frequencies
# number of upper harmonics to consider
psd_n_harms = 20
# number of neighbouring frequency bins to consider per side for SNR
# baseline correction with psd_base_bins will be applied for z-score
psd_snr_bins = 10
# number of neighbouring frequency bins (per side)
# for "baseline-correction" of PSDs (can be 0)
psd_base_bins = 10
# number of bins as "gap" between neighours (n_bins) and target frequency
psd_n_gap = 0

###############################################################################
# Maxfilter etc.

# which files to maxfilter and how to name them after sss
# [before maxfilter], [after maxfilter], [condition labels],
# [presentation/oddball frequencies]

sss_map_fnames = {
    1: (['rest1_raw', 'rest2_raw',
         'Faces_raw',
         'hfwpw1_raw', 'hfwpw2_raw', 'HFWPW3_raw',
         'LFWPW1_raw', 'LFWPW2_raw', 'LFWPW3_raw',
         'HFWLFW1_raw', 'HFWLFW2_raw', 'HFWLFW3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhfw1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    2: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    3: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    4: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
}


# parameters for Neuromag maxfilter command
MF = {
     'NM_cmd': '/imaging/local/software/neuromag/bin/util/maxfilter-2.2.12',
     'cal': '/neuro/databases/sss/sss_cal.dat',
     'ctc': '/neuro/databases/ctc/ct_sparse.fif',
     'st_duration': 10.,
     'st_correlation': 0.98,
     'origin': (0., 0., 0.045),
     'in': 8,
     'out': 3,
     'regularize': 'in',
     'frame': 'head',
     'mv': 'inter',
     'trans': 6}  # which file to use for -trans

# for correcting EEG electrode positions
check_cmd = '/imaging/local/software/mne/mne_2.7.3/x86_64/\
MNE-2.7.3-3268-Linux-x86_64//bin/mne_check_eeg_locations \
--file %s --fix'

### FILTERING, EVENTS

# define the stim channel
stim_channel = 'STI101'

# bandpass filter frequencies
l_freq, h_freq = None, 40.

# occipital EEG electrodes
occ_eeg = ['EEG0%.2d' % i for i in list(range(40, 61)) +
           list(range(65, 75))]

###
# TFR
###

tfr = {}
tfr['epoch'] = [-3., 60.]  # latency range of epoch for TFR (s)
tfr['freqs'] = np.arange(0.5, 15.25, 0.25)  # frequencies for TFR (Hz)

########################################################
# Edited for FPVS up to here
########################################################

### Epoching, Averaging

# stimulus projector delay
delay = 0.0345

# separate triggers for target detection and localiser tasks
event_id = {}

event_id['td'] = {
        'stim': 254, # beginning of grating/sound
        'aud': 101, # trigger sent on presentation of the word
        'vis': 102,
        'act': 103,
        'abs': 104
        # 'catch': 105,
        # 'filler': 106
    }

event_id['loc'] = {'vis': 107, 'aud': 108}

# event_id['av'] = {'stim': 254}

# epoch interval and baseline
epo_t = {'loc': {'epo': [-0.5, 1.2], 'baseline': (-0.5, 0.)},
         # 'td': {'epo': [-1.2, 0.5], 'baseline': (-1.2, -0.7)},
         'td': {'epo': [-1.2, 0.5], 'baseline': (-0.3, 0.0)},
         'av': {'epo': [-0.5, 1.2], 'baseline': (-0.5, 0.)}
         }

# mark which subjects have an issue with the ICA --> don't plot the EOG epochs

bad_ica = [5, 13, 16, 18, 23, 24]


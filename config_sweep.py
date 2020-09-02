"""
===========
Config file for FPVS with Frequency Sweep
===========
"""

import os
from os import path as op

import sys

import numpy as np

###############################################################################

# IDs of subjects to process (SLURM and Grand-Average)
do_subjs = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18]
# do_subjs = [18]
# removed:
# 6: no MRI
# 15: MEG artifacts due to titanium plate
# 17: no MRI

# which conditions to process (after Maxfilter and filtering)
do_conds = ['face']

# paths to data:

# path to acquired raw data
cbu_path = '/megdata/cbu/fpvs'

# path to data for pre-processing
data_path = '/group/erp/data/olaf.hauk/MEG/FPVS/data_Federica'

# path to Freesurfer-preprocessed MRIs
subjects_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data_Federica/MRI'

# for grand-mean results
grandmean_path = '/group/erp/data/olaf.hauk/MEG/FPVS/data_Federica/GM'

# for data exported to ASCII, Matlab etc.
export_path = '/group/erp/data/olaf.hauk/MEG/FPVS/data_Federica/export'
do_export = True  # whether to export raw spectra to Matlab or not

if not os.path.isdir(data_path):  # create if necessary
    os.mkdir(data_path)

# Compute TRF (1) or not (0)
do_tfr = 0

# use ICAed files or not; end of filenames
raw_ICA_suff = 'ica_raw'
# raw_ICA_suff = 'raw'

###############################################################################
# Mapping betwen filenames and subjects

map_subjects = {
    1: ('meg19_0380', '191008'),  # pilot frequency sweep
    2: ('meg19_0381', '191015'),  
    3: ('meg19_0383', '191018'),  
    4: ('meg19_0384', '191021'),
    5: ('meg19_0389', '191022'),    
    6: ('meg19_0391', '191028'),
    7: ('meg19_0392', '191029'),  
    8: ('meg19_0396', '191031'),  # ARTEFACTS/BAD MEG CHANNELS due to zip on the participant's top 
    9: ('meg19_0400', '191105'),  # in faces: pressing button with left hand (noticed in 2nd trial)
    10: ('meg19_0406', '191108'),   #slouched down several times during the recording
    11: ('meg19_0407', '191111'),   
    12: ('meg19_0412', '191114'),   
    13: ('meg19_0414', '191115'),   #big head + sleepy
    14: ('meg19_0417', '191118'),   #sleepy
    15: ('meg19_0421', '191122'),   # Titanium plate, ***************check if need to DISCARD noisy after Maxfilter probably discard
    16: ('meg19_0422', '191122'),
    17: ('meg19_0425', '191125'),   #sleepy but NO alpha + some strange blinks
    18: ('meg19_0442', '191205')    
}

# subject names of MRI data
mri_subjects = {
    1: ('CBU000000'),
    2: ('CBU180712'),
    3: ('CBU190847'),
    4: ('CBU190899'),
    5: ('CBU200028'),
    6: (''),  # didn't respond
    7: ('CBU190829'),
    8: ('CBU190940'),
    9: ('CBU190917'),
    10: ('CBU190200'),
    11: ('CBU190893'),
    12: ('CBU190613'),
    13: ('CBU190840'),
    14: ('CBU170707'),
    15: ('CBU200023'),
    16: ('CBU190994'),
    17: (''),  # to be recorded
    18: ('CBU200027')
}


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
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
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
    5: (['Rest1_raw', 'Rest2_avg',
         'Faces_raw',
         'PWHF1_raw', 'PWHF2_raw', 'PWHF3_raw',
         'PWLF1_raw', 'PWLF2_raw', 'PWLF3_raw',
         'LFHF1_raw', 'LFHF2_raw', 'LFHF3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    6: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw1', 'lfhf2_raw', 'lfhf3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    7: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    8: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf21_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    9: (['rest1_raw', 'rest2_raw',
         'Faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raew', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
        ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    10: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    11: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_araw', 'lfhf2_raw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    12: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    13: (['rest_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_avg', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    14: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    15: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    16: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    17: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_raw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw']),
    18: (['rest1_raw', 'rest2_raw',
         'faces_raw',
         'pwhf1_raw', 'pwhf2_raw', 'pwhf3_raw',
         'pwlf1_raw', 'pwlf2_raw', 'pwlf3_raw',
         'lfhf1_raw', 'lfhf2_rw', 'lfhf3_raw'],
         ['rest1_sss_raw', 'rest2_sss_raw',
         'faces_sss_raw',
         'pwhf1_sss_raw', 'pwhf2_sss_raw', 'pwhf3_sss_raw',
         'pwlf1_sss_raw', 'pwlf2_sss_raw', 'pwlf3_sss_raw',
         'lfhf1_sss_raw', 'lfhf2_sss_raw', 'lfhf3_sss_raw'])
}


###############################################################################
# Bad channels

bad_channels = {
    1: {'eeg': ['EEG028'],
        'meg': ['MEG1123', 'MEG2223', 'MEG0813']},
    2: {'eeg': ['EEG041'],
        'meg': ['MEG1123', 'MEG0813', 'MEG1412']},
    3: {'eeg': [],
        'meg': ['MEG2312', 'MEG2311', 'MEG0413']},
    4: {'eeg': ['EEG012', 'EEG074'],
        'meg': ['MEG0723', 'MEG1123', 'MEG0813']},
    5: {'eeg': ['EEG029', 'EEG039', 'EEG050', 'EEG056', 'EEG071'],  # also 'EEG049' is a bit noisy
        'meg': ['MEG1711', 'MEG2312', 'MEG0813', 'MEG2311']},
    6: {'eeg': ['EEG023', 'EEG040', 'EEG043'],
        'meg': ['MEG2323', 'MEG0813']},
    7: {'eeg': ['EEG019','EEG023'],
        'meg': ['MEG0233', 'MEG1121', 'MEG0813', 'MEG1123']},
    8: {'eeg': ['EEG004', 'EEG018', 'EEG029', 'EEG039', 'EEG050'],
        'meg': []},
    9: {'eeg': [],
        'meg': ['MEG0813']},
    10: {'eeg': [],
        'meg': ['MEG1131']},
    11: {'eeg': [],
        'meg': ['MEG0813']},
    12: {'eeg': ['EEG029', 'EEG048'],
        'meg': ['MEG0813']},
    13: {'eeg': ['EEG039', 'EEG050'],
        'meg': ['MEG1222']},
    14: {'eeg': ['EEG032', 'EEG045', 'EEG047'],           # need to check how it looks with this channel interpolation
        'meg': ['MEG2642', 'MEG0412', 'MEG1713', 'MEG2312', 'MEG2323']},
    15: {'eeg': [],
        'meg': ['MEG2511', 'MEG0813', 'MEG0933']},
    16: {'eeg': ['EEG045'],
        'meg': ['MEG1322', 'MEG2223']},
    17: {'eeg': ['EEG039', 'EEG048', 'EEG050', 'EEG055'],
        'meg': ['MEG0813', 'MEG1712']},
    18: {'eeg': ['EEG008', 'EEG021', 'EEG029', 'EEG045'],
        'meg': ['MEG2323', 'MEG1943', 'MEG1741']},

}


# create subject-specific data directories if necessary
for ss in map_subjects:
    # subject-specific sub-dir, e.g. maxfiltered raw data
    subj_dir = op.join(data_path, map_subjects[ss][0])
    if not op.isdir(subj_dir):
        print('Creating directory %s.' % subj_dir)
        os.mkdir(subj_dir)
    # subject-specific sub-dir for evoked data
    subj_dir_ave = op.join(data_path, map_subjects[ss][0], 'AVE')
    if not op.isdir(subj_dir_ave):
        print('Creating directory %s.' % subj_dir_ave)
        os.mkdir(subj_dir_ave)
    # subject-specific sub-dir for epochs
    subj_dir_epo = op.join(data_path, map_subjects[ss][0], 'EPO')
    if not op.isdir(subj_dir_epo):
        print('Creating directory %s.' % subj_dir_epo)
        os.mkdir(subj_dir_epo)
    # subject-specific sub-dir for source space data
    subj_dir_stc = op.join(data_path, map_subjects[ss][0], 'STC')
    if not op.isdir(subj_dir_stc):
        print('Creating directory %s.' % subj_dir_stc)
        os.mkdir(subj_dir_stc)
    fig_dir = op.join(data_path, map_subjects[ss][0], 'Figures')  # subject figure dir
    if not op.isdir(fig_dir):
        print('Creating directory %s.' % fig_dir)
        os.mkdir(fig_dir)
    fig_dir = op.join(data_path, map_subjects[ss][0], 'Figures_ICA')  # subject figure dir
    if not op.isdir(fig_dir):
        print('Creating directory %s.' % fig_dir)
        os.mkdir(fig_dir)

if not op.isdir(grandmean_path):
    os.mkdir(grandmean_path)
    os.mkdir(op.join(grandmean_path, 'Figures'))
    os.mkdir(op.join(grandmean_path, 'Figures_ICA'))
if not op.isdir(op.join(grandmean_path, 'AVE')):
    os.mkdir(op.join(grandmean_path, 'AVE'))
    os.mkdir(op.join(grandmean_path, 'STC'))

# For subjects without clean ECG channel,
# use the following magnetometers in ICA (else specify '' to use ECG)
ECG_channels = {
    1: '',
    2: '',
    3: '',
    4: '',
    5: '',
    6: '',
    7: '',
    8: '',
    9: '',
    10: '',
    11: '',
    12: '',
    13: '',
    14: '',
    15: '',
    16: '',
    17: '',
    18: ''
}

# Artefact rejection thresholds
# for ICA, covariance matrix
reject = dict(grad=4e-10, mag=1e-11, eeg=1e-3)

###############################################################################
# ERPs

# artefact rejection thresholds for epoching
epo_reject = dict(grad=4e-10, mag=1e-11, eeg=1e-3)

# baseline in s
epo_baseline = (-.2, 0.)

# epoch interval in s
epo_t1, epo_t2 = -.2, .5

###############################################################################
# FPVS

# resample to this sampling frequency before computing PSDs
# only resample if not None
psd_resample = None

# window size of FFT
# determined with mne.filter.next_fast_len()
# chosen to match one sweep duration in the word conditions
psd_nfft = 12000

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
# depends on frequency range used for PSD
fpvs_n_harms_odd = 10
fpvs_n_harms_base = 10

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
topo_times['words'] = np.array([1., 2., 3., 4., 6., 10., 12.])
topo_times['faces'] = np.array([1.2, 2.4, 3.6, 4.8, 6., 12.])

# Frequency range for PSD (will determine number of possible harmonics)
psd_fmin = 0.
psd_fmax = 140.

# for plot_joint() display
# These are actually frequencies in Hz, for PSD plots
crop_times = [0.1, 30.]

# number of neighbouring frequency bins to consider per side for SNR
# baseline correction with psd_snr_bins will be applied for z-score
# depends on frequency resolution
psd_snr_bins = {}
psd_snr_bins['faces'] = 8  # about 0.8 Hz
psd_snr_bins['words'] = 8  # about 0.8 Hz
# number of neighbouring frequency bins (per side)
# number of bins as "gap" between neighours (n_bins) and target frequency
psd_n_gap = 1
# number of peak channels to select for plots
n_peak = 4

###############################################################################
# Maxfilter etc.

# parameters for Neuromag maxfilter command
# Make sure to use Vectorview files!
MF = {
     'NM_cmd': '/imaging/local/software/neuromag/bin/util/maxfilter-2.2.12',
     'cal': '/neuro/databases_vectorview/sss/sss_cal.dat',
     'ctc': '/neuro/databases_vectorview/ctc/ct_sparse.fif',
     'st_duration': 10.,
     'st_correlation': 0.98,
     'origin': (0., 0., 0.045),
     'in': 8,
     'out': 3,
     'regularize': 'in',
     'frame': 'head',
     'mv': 'inter',
     'trans': 6}  # which file to use for -trans within subject

# for correcting EEG electrode positions
check_cmd = '/imaging/local/software/mne/mne_2.7.3/x86_64/\
MNE-2.7.3-3268-Linux-x86_64//bin/mne_check_eeg_locations \
--file %s --fix'

### FILTERING, EVENTS

# define the stim channel
stim_channel = 'STI101'

# bandpass filter frequencies
l_freq, h_freq = 0.1, 140.

# from Retter et al., bioRxiv 2019, p. 10f.:
# "OT": right: channels P10; P8; PO8; PO10; PO12;
#       and left: P9; P7; PO7; PO9; PO11
# "MO": O2; POI2; I2; Iz; OIz; Oz; POOz; O1; POI1; I1
electrode_ROIs = {}
electrode_ROIs['OT_R'] = ['EEG065', 'EEG060', 'EEG070', 'EEG003']  # we don't have PO12
electrode_ROIs['OT_L'] = ['EEG051', 'EEG052', 'EEG066', 'EEG001']
electrode_ROIs['MO'] = ['EEG073', 'EEG074', 'EEG072', 'EEG071']  # we only have O2, IZ, Oz, O1

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

### Source Space
stc_morph = 'fsaverage'

# vertex size
src_spacing = 5
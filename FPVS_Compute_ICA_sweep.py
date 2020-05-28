#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Compute ICA for FPVS Frequency Sweep.

Based on Fiff_Compute_ICA.py.
==========================================
OH July 2019
"""

# TO DO: change to compute ICA across concatenated raw files/epochs
# often only ~6 blinks per run
# Subject 18 has not detectable eye-blinks in at least one run, lfhf1

import sys

from os import path as op
import numpy as np

# import matplotlib
# matplotlib.use('Agg') #  for running graphics on cluster ### EDIT
from matplotlib import pyplot as plt

from importlib import reload

import mne
from mne.preprocessing import ICA, create_eog_epochs, create_ecg_epochs
from mne.report import Report

import config_sweep as config
reload(config)

print('MNE Version: %s\n\n' % mne.__version__)  # just in case
print(mne)

# whether to show figures on screen or just write to file
show = False

# "emulate" the args from ArgParser in Fiff_Compute_ICA.py
# filenames depend on subject, the rest are variables
class create_args:
  def __init__(self, FileRaw, FileICA, FileHTML):
    self.FileRaw = FileRaw
    self.FileICA = FileICA
    self.FileHTML = FileHTML
    self.EOG = ['EOG062']
    self.ECG = ['ECG063']  # test 'synth' option if no ECG channel present
    self.maxEOG = 2
    self.maxECG = 2

    self.ECGmeth = 'ctps'
    self.EOGthresh = 2.5
    self.ECGthresh = 0.05

    self.ChanTypes = ['eeg', 'meg']
    self.RejEEG = config.reject['eeg']
    self.RejGrad = config.reject['grad']
    self.RejMag = config.reject['mag']
    self.n_pca_comps = '0.99'  # string required
    self.method = 'infomax'


def run_Compute_ICA(sbj_id):
    """Compute ICA for one subject."""
    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    # Concatenate raws. These raw files are not huge.
    raw_all = []  # will contain list of all raw files
    for raw_stem_in in sss_map_fname[1]:

        FileRaw = op.join(sbj_path, raw_stem_in[:-7] + 'sss_f_raw')

        # Read raw data info
        raw1 = mne.io.read_raw_fif(FileRaw + '.fif', preload=True)

        raw_all.append(raw1)

    # concatenate raws
    print('Concatenating %d raw files.' % len(raw_all))
    raw = mne.concatenate_raws(raw_all)

    # -ica.fif will be appended
    FileICA = op.join(sbj_path, config.map_subjects[sbj_id][0] + '_sss_f_raw')

    # -ica.html will be appended
    FileHTML = op.join(sbj_path, config.map_subjects[sbj_id][0] + '_sss_f_raw')

    # define variables for the following ICA pipeline
    # this would be from command line arguments of Fiff_Compute_ICA.py
    args = create_args(FileRaw, FileICA, FileHTML)

    # If a channel for ECG detection explicity specified, use it
    if config.ECG_channels[sbj_id] != '':

        args.ECG = [config.ECG_channels[sbj_id]]

    # otherwise use ECG from data, but if not present, dont' do ICA for ECG
    elif not raw.__contains__('ecg'):

        args.ECG = []

        print('###\nNo ECG found in raw data, so I am not doing it.\n###')

    # Now turn the "fake command line parameters" into variables for the
    # analysis pipeline

    # if float, select n_components by explained variance of PCA
    if '.' in args.n_pca_comps:

        n_components = float(args.n_pca_comps)
        print('Number of PCA components by fraction of variance (%f)' %
              n_components)

    else:

        n_components = int(args.n_pca_comps)
        print('Number of PCA components: %d.' % n_components)

    method = args.method  # for comparison with EEGLAB "extended-infomax"
    print('\nUsing ICA method %s.' % method)

    decim = 3  # downsample data to save time

    # same random state for each ICA (not sure if beneficial?)
    random_state = 23

    # raw data input filename, not needed here for concatenated raws
    raw_fname_in = args.FileRaw + '.fif'

    # filename for ICA output
    if args.FileICA == '':

        ica_fname_out = args.FileRaw + '-ica.fif'

    else:

        ica_fname_out = args.FileICA + '-ica.fif'


    # filename for ICA output
    if args.FileHTML == '':

        fname_html = args.FileRaw + '-ica.html'

    else:

        fname_html = args.FileHTML + '-ica.html'

    ###
    # START ICA
    ###

    report = Report(subject=config.map_subjects[sbj_id][0], title='ICA:')

    # print('###\nReading raw file %s.' % raw_fname_in)

    # # Read raw data
    # raw = mne.io.read_raw_fif(raw_fname_in, preload=True)

    # check if EEG in raw data
    if not raw.__contains__('eeg'):

        args.ChanTypes = ['meg']

        print('###\nNo EEG found in raw data, continuing with MEG only.\n###')

    # They say high-pass filtering helps
    print('Band-pass filtering raw data between 1 and 40 Hz.')
    raw.filter(1., 40., fir_design='firwin')

    # which channel types to use
    to_pick = {'meg': False, 'eeg': False, 'eog': False, 'stim': False,
               'exclude': 'bads'}

    # pick channel types as specified
    print('Using channel types: ')
    for chtype in args.ChanTypes:

        print(chtype + ' ')
        to_pick[chtype.lower()] = True

    picks_meg_eeg_eog = mne.pick_types(raw.info, meg=to_pick['meg'],
                                       eeg=to_pick['eeg'],
                                       eog=True, ecg=True,
                                       stim=to_pick['stim'],
                                       exclude=to_pick['exclude'])

    # to remove non-physiological artefacts (parameters based on MNE example)
    reject = {}
    if to_pick['meg'] == True:

        reject['mag'] = args.RejMag
        reject['grad'] = args.RejGrad
        print('Thresholds for MEG: Grad %.1e, Mag %.1e.' % (reject['grad'], reject['mag']))

    if to_pick['eeg'] == True:

        reject['eeg'] = args.RejEEG
        print('Threshold for EEG: %.1e.' % reject['eeg'])

    picks_meg = mne.pick_types(raw.info, meg=to_pick['meg'],
                               eeg=to_pick['eeg'],
                               eog=to_pick['eog'], stim=to_pick['stim'],
                               exclude=to_pick['exclude'])

    # Compute ICA model ###################################################

    print('###\nDefine the ICA object instance using %s. \
          Number of PCA components based on: %s.' %
          (method, str(n_components)))

    ica = ICA(n_components=n_components, method=method,
              random_state=random_state)

    print('Fitting ICA.')

    ica.fit(raw, picks=picks_meg, decim=decim, reject=reject)
    print(ica)

    print('Plotting ICA components.')

    # plot for specified channel types
    for ch_type in reject.keys():

        fig_ic = ica.plot_components(ch_type=ch_type, show=show)

        captions = [ch_type.upper() + ' Components' for i in fig_ic]

        report.add_figs_to_section(fig_ic, captions=captions,
                                   section='ICA Components', scale=1)

    # indices of ICA components to be removed across EOG and ECG
    ica_inds = []

    ###
    # EOG COMPONENTS
    ###

    # for all specified EOG channels
    eog_inds = []  # ICA components found to be bad for EOG
    eog_scores = []  # corresponding ICA scores

    for eog_ch in args.EOG:

        print('\n###\nFinding components for EOG channel %s.\n' % eog_ch)

        # get single EOG trials
        eog_epochs = create_eog_epochs(raw, ch_name=eog_ch, reject=reject)

        eog_average = eog_epochs.average()  # average EOG epochs

        # find via correlation
        inds, scores = ica.find_bads_eog(eog_epochs, ch_name=eog_ch,
                                         threshold=args.EOGthresh)

        if inds != []:  # if some bad components found

            print('###\nEOG components and scores for channel %s:\n' % eog_ch)
            for [ee, ss] in zip(inds, scores):

                print('%d: %.2f\n' % (ee, ss))

            # look at r scores of components
            fig_sc = ica.plot_scores(scores, exclude=inds, show=show)

            report.add_figs_to_section(fig_sc, captions='%s Scores' %
                                       eog_ch, section='%s ICA component \
                                       scores' % eog_ch, scale=1)

            print('Plotting raw ICA sources.')
            fig_rc = ica.plot_sources(raw, show=show)

            report.add_figs_to_section(fig_rc, captions='%s Sources' %
                                       eog_ch, section='%s raw ICA sources'
                                       % eog_ch, scale=1)

            print('Plotting EOG average sources.')
            # look at source time course
            fig_so = ica.plot_sources(eog_average, show=show)

            report.add_figs_to_section(fig_so, captions='%s Sources' %
                                       eog_ch, section='%s ICA Sources' %
                                       eog_ch, scale=1)

            print('Plotting EOG epochs properties.')
            fig_pr = ica.plot_properties(eog_epochs, picks=inds,
                                         psd_args={'fmax': 35.},
                                         image_args={'sigma': 1.},
                                         show=show)

            txt_str = '%s Properties' % eog_ch
            captions = [txt_str for i in fig_pr]

            report.add_figs_to_section(fig_pr, captions=captions,
                                       section='%s ICA Properties' %
                                       eog_ch, scale=1)

            print(ica.labels_)

            # Remove ICA components #######################################
            fig_ov = ica.plot_overlay(eog_average, exclude=inds, show=show)
            # red -> before, black -> after.

            report.add_figs_to_section(fig_ov, captions='%s Overlay' %
                                       eog_ch, section='%s ICA Overlay' %
                                       eog_ch, scale=1)

            plt.close('all')

            eog_inds += inds  # keep bad ICA components

            # keep scores for bad ICA components
            eog_scores += list(scores[inds])

        else:

            print('\n!!!Nothing bad found for %s!!!\n' % eog_ch)

    if eog_inds != []:  # if there are bad ECG components

        # deal with case where there are more bad ICA components than
        # specified
        n_comps = np.min([args.maxEOG, len(eog_inds)])

        print('\n###\nUsing %d out of %d detected ICA components for EOG.'
              % (n_comps, len(eog_inds)))

        for [c, s] in zip(eog_inds, eog_scores):

            print('Component %d with score %f.' % (c, s))
        # sort to find ICA components with highest scores
        idx_sort = np.argsort(np.abs(eog_scores))

        # only keep desired number of bad ICA components with high scores
        ica_inds += [eog_inds[idx] for idx in idx_sort[-n_comps:]]

    #
    # ECG COMPONENTS
    #

    # for all specified EOG channels

    ecg_inds = []  # ICA components found to be bad for ECG
    ecg_scores = []  # corresponding ICA scores

    for ecg_ch in args.ECG:

        if ecg_ch == 'synth':

            print('Creating synthetic ECG channel.')

            # check which channel, if any, is ECG
            ecg_idx = np.where(['ECG' in ch for ch in
                               raw.info['ch_names']])[0]

            # if there is an ECG channel, change it
            if not ecg_idx.shape[0] == 0:

                ecg_name = raw.info['ch_names'][ecg_idx[0]]

                raw.set_channel_types({ecg_name: 'misc'})

            # create synthetic ECG channel across MEG channels
            ecg_ch_name = None  # for create_ecg_epochs
            ecg_find_name = 'ECG-SYN'
            keep_ecg = True

        else:

            ecg_ch_name = ecg_ch
            ecg_find_name = ecg_ch
            keep_ecg = False

        print('\n###\nFinding components for ECG channel %s.\n' % ecg_ch)

        # get single ECG trials
        ecg_epochs = create_ecg_epochs(raw, ch_name=ecg_ch_name,
                                       keep_ecg=keep_ecg, reject=reject)

        ecg_average = ecg_epochs.average()  # average ECG epochs

        # find via cross-trial phase statistics
        inds, scores = ica.find_bads_ecg(ecg_epochs, ch_name=ecg_find_name,
                                         method=args.ECGmeth,
                                         threshold=args.ECGthresh)

        if inds != []:  # if some bad components found

            print('ECG components and scores:\n')
            for [ee, ss] in zip(inds, scores):

                print('%d: %.2f\n' % (ee, ss))

            # look at r scores of components
            fig_sc = ica.plot_scores(scores, show=show)

            report.add_figs_to_section(fig_sc, captions='%s Scores' %
                                       ecg_ch, section='%s component \
                                       scores' % ecg_ch, scale=1)

            print('Plotting raw ICA sources.')
            fig_rc = ica.plot_sources(raw, show=show)

            report.add_figs_to_section(fig_rc, captions='%s Sources' %
                                       ecg_ch, section='%s raw sources'
                                       % ecg_ch, scale=1)

            print('Plotting ECG average sources.')
            # look at source time course
            fig_so = ica.plot_sources(ecg_average, show=show)

            report.add_figs_to_section(fig_so, captions='%s Sources' %
                                       ecg_ch, section='%s ICA Sources' %
                                       ecg_ch, scale=1)

            print('Plotting ECG epochs properties.')
            fig_pr = ica.plot_properties(ecg_epochs, picks=inds,
                                         psd_args={'fmax': 35.},
                                         image_args={'sigma': 1.},
                                         show=show)

            txt_str = '%s Properties' % ecg_ch
            captions = [txt_str for i in fig_pr]

            report.add_figs_to_section(fig_pr, captions=captions,
                                       section='%s ICA Properties' %
                                       ecg_ch, scale=1)

            print(ica.labels_)

            # Remove ICA components #######################################
            fig_ov = ica.plot_overlay(ecg_average, exclude=inds, show=show)
            # red -> before, black -> after. Yes! We remove quite a lot!

            report.add_figs_to_section(fig_ov, captions='%s Overlay' %
                                       ecg_ch, section='%s ICA Overlay' %
                                       ecg_ch, scale=1)

            plt.close('all')

            ecg_inds += inds  # keep bad ICA components
            ecg_scores += list(scores[inds])  # keep bad ICA components

        else:

            print('\n!!!Nothing bad found for %s!!!\n' % ecg_ch)

    if ecg_inds != []:  # if there are bad ECG components

        # deal with case where there are more bad ICA components than
        # specified
        n_comps = np.min([args.maxECG, len(ecg_inds)])

        print('\n###\nUsing %d out of %d detected ICA components for ECG.'
              % (n_comps, len(ecg_inds)))

        for [c, s] in zip(ecg_inds, ecg_scores):

            print('Component %d with score %f.' % (c, s))
        # sort to find ICA components with highest scores
        idx_sort = np.argsort(np.abs(ecg_scores))

        # only keep desired number of bad ICA components with high scores
        ica_inds += [ecg_inds[idx] for idx in idx_sort[-n_comps:]]

    if ica_inds == []:

        print('###\nNo bad components found anywhere.\n###')

    # specify components to be removed
    ica.exclude = ica_inds

    ###
    # SAVE ICA
    ###

    # from now on the ICA will reject this component even if no exclude
    # parameter is passed, and this information will be stored to disk
    # on saving

    print('\nSaving ICA to %s' % (ica_fname_out))
    ica.save(ica_fname_out)

    print('Saving HTML report to {0}'.format(fname_html))
    report.save(fname_html, overwrite=True, open_browser=True)


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0,len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    run_Compute_ICA(ss)

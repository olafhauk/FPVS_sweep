#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Compute grand-average of the outputs of FPVS_PSD_sweep.py.

In order to average across all subjects, don't specify arguments, or specify
a number larger than the largest subject ID (e.g. using using SLURM).
For example:
run FPVS_GrandAverage_PSDs
run FPVS_GrandAverage_PSDs 1 2 3 4
run FPVS_GrandAverage_PSDs 99
==========================================

OH, February 2020
"""

# To do:
# data are not in format for easy averaging across subjects
# then plotting needs changes accordingly

# for Evoked data are in one file for all frequencies
# for STC data are in separate files per condition and freq

import sys

from os import path as op

import numpy as np

from scipy.stats import ttest_rel

from matplotlib import pyplot as plt

from copy import deepcopy

from importlib import reload

import mne
from mne.report import Report
from mne.source_estimate import SourceEstimate
from mne.evoked import EvokedArray

import FPVS_functions
reload(FPVS_functions)

from FPVS_functions import grand_average_evoked_arrays

import config_sweep as config
reload(config)

import FPVS_functions as Ff
reload(Ff)


print(mne.__version__)

# perform TFR of raw data or not
# do_tfr = config.do_tfr

# sub-directory for figures per subject
# separate for ICAed and non-ICAed data
if 'ica' in config.raw_ICA_suff:
    figs_dir = 'Figures_ICA'
else:
    figs_dir = 'Figures'

close_fig = 1  # close figures only if close_fig==1

# plt.ion() # interactive plotting

# for some plots of SNRs
unit_scalings = dict(eeg=1., mag=1., grad=1.)

# Base frequencies for frequency sweep for words (not faces)
freqs_all = [str(ff) for ff in config.fpvs_freqs]

# separate filename prefixes for ICAed and non-ICAed data
prefix = ''
if 'ica' in config.raw_ICA_suff:
    prefix = 'ICA'

# Which modalities and results to process

# all psd results for evoked and STC
# individual subjects and GM
modals = ['evo', 'stc']
gm_modals = ['evo_gm', 'stc_gm']
# modals = ['stc']
# gm_modals = ['stc_gm']

# for evoked
types = ['psd', 'psd_z', 'psd_sum_odd', 'psd_sum_base', 'psd_harm_odd',
         'psd_harm_base', 'psd_harm_topos_odd', 'psd_harm_topos_base']

# only for evoked: data for peak channels per condition
evo_types = [
    'peak_odd', 'z_peak_odd', 'harm_odd_peak_odd',
    'harm_base_peak_odd', 'peak_base', 'z_peak_base', 'harm_odd_peak_base',
    'harm_base_peak_base', 'peak_harm_topos_odd', 'peak_harm_topos_base']

# for STCs
stc_types = ['psd', 'psd_sum_odd', 'psd_sum_base', 'psd_harm_odd',
             'psd_harm_base', 'psd_harm_topos_odd', 'psd_harm_topos_base']

# conditions
# conds = ['face', 'pwhf', 'pwlf', 'lfhf']
conds = config.do_conds


def grand_average_psds(sbj_ids):
    """Grand-average PSDs and derivatives across sbj_ids."""
    # initialise html report for one subject

    print('Grand-averaging subjects:')
    print(*sbj_ids)

    # report = Report(subject='GM', title='FPVS PSDs GM')

    # # get condition names and frequency names
    # # assumed to be consistent across subjects
    # sss_map_fname = config.sss_map_fnames[sbj_ids[0]]
    # conds = []  # names of conditions
    # for raw_stem_in in sss_map_fname[1][2:]:

    #     conds.append(raw_stem_in[:4])

    # conds = np.unique(conds)

    # initialise

    psds = {}  # individual subjects and GM

    do_modals = modals + gm_modals

    for modal in do_modals:

        print(modal)

        psds[modal] = {}  # type of data

        do_types = types

        if modal == 'evo':  # add other types

            do_types = do_types + evo_types

        for tt in do_types:

            psds[modal][tt] = {}  # type of processed PSD

            for cond in conds:

                psds[modal][tt][cond] = {}  # sweep frequencies

                if cond == 'face':  # no frequency sweep for faces

                    freqs = ['6.0']  # base frequency for this condition (Hz as string)

                else:  # for all word condition, use all sweep frequencies

                    # base frequencies for this condition (Hz as string)
                    freqs = freqs_all

                for freq in freqs:

                    psds[modal][tt][cond][freq] = []  # subjects

    # initialise array for electrode ROIs for group statistics
    roi_elecs = {}
    for roi in config.electrode_ROIs:

        roi_elecs[roi] = {}

        for cond in conds:

            roi_elecs[roi][cond] = {}

            if cond == 'face':  # no frequency sweep for faces

                freqs = ['6.0']  # base frequency for this condition (Hz as string)

            else:  # for all word condition, use all sweep frequencies

                # base frequencies for this condition (Hz as string)
                freqs = freqs_all

            for freq in freqs:

                roi_elecs[roi][cond][freq] = {'odd': np.zeros(len(sbj_ids)),
                                              'base': np.zeros(len(sbj_ids))}

    # Reading evoked data, getting data for channel groups
    if 'evo' in modals:

        print('Reading evoked data.')

        modal = 'evo'

        for cond in conds:  # conditions

            print('Condition: %s.' % cond)

            if cond == 'face':  # no frequency sweep for faces

                freqs = ['6.0']  # base frequency for condition (Hz as str)

                freq_odd = 1.2  # oddball frequency for this condition (Hz)

            else:  # for all word condition, use all sweep frequencies

                # base frequencies for this condition (Hz as string)
                freqs = freqs_all

                freq_odd = 1.0  # oddball frequency the same for all sweeps

            # for Evoked data are in one file for all frequencies
            # for STC data are in separate files per condition and freq
            for [ss, sbj_id] in enumerate(sbj_ids):  # across all subjects

                # path to subject's data
                sbj_path = op.join(config.data_path,
                                   config.map_subjects[sbj_id][0])

                print('Reading PSD results from evoked files:')

                # Read Evoked

                fname_evo =\
                    op.join(sbj_path, 'AVE', '%sPSDTopo_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd = mne.read_evokeds(fname_evo)

                fname_evo =\
                    op.join(sbj_path, 'AVE', '%sPSDTopoZ_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd_z = mne.read_evokeds(fname_evo)

                fname_evo =\
                    op.join(sbj_path, 'AVE', '%sPSDHarm_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd_harm_odd = mne.read_evokeds(fname_evo)

                fname_evo = \
                    op.join(sbj_path, 'AVE', '%sPSDHarmBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd_harm_base = mne.read_evokeds(fname_evo)

                fname_evo =\
                    op.join(sbj_path, 'AVE', '%sPSDSumTopoOdd_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd_sum_odd = mne.read_evokeds(fname_evo)

                fname_evo =\
                    op.join(sbj_path, 'AVE', '%sPSDSumTopoBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd_sum_base = mne.read_evokeds(fname_evo)

                fname_evo =\
                    op.join(sbj_path, 'AVE', '%sPSDSumToposOdd_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd_harm_topos_odd = mne.read_evokeds(fname_evo)

                fname_evo =\
                    op.join(sbj_path, 'AVE', '%sPSDSumToposBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
                print(fname_evo)
                psd_harm_topos_base = mne.read_evokeds(fname_evo)

                for (fi, freq) in enumerate(freqs):

                    print('Frequency: %s.' % freq)

                    psds[modal]['psd'][cond][freq].append(psd[fi])

                    psds[modal]['psd_z'][cond][freq].append(psd_z[fi])

                    psds[modal]['psd_sum_odd'][cond][freq].append(psd_sum_odd[fi])

                    psds[modal]['psd_sum_base'][cond][freq].append(psd_sum_base[fi])

                    psds[modal]['psd_harm_topos_odd'][cond][freq].append(psd_harm_topos_odd[fi])

                    psds[modal]['psd_harm_topos_base'][cond][freq].append(psd_harm_topos_base[fi])

                    psds[modal]['psd_harm_odd'][cond][freq].append(psd_harm_odd[fi])

                    psds[modal]['psd_harm_base'][cond][freq].append(psd_harm_base[fi])

                    # hack, float-to-string-to-float-again
                    # to be consistent with FPVS_PSD_sweep_plot.py
                    basefreq = float(freq)

                    # Get max channels from z-scored PSD at base frequency
                    # not oddball frequency, which would be biased.
                    # This evoked is for condition cond, current subject and
                    # current frequency freq.
                    evoked = deepcopy(psd_z[fi])

                    # Find channels with maximum Z-scores per channel type
                    # for base frequency
                    # "Latency" is frequency in Hz divided by 1000
                    peak_times_base = [basefreq]

                    peak_ch_types_base = Ff.peak_channels_evoked(
                        evoked=evoked, peak_times=peak_times_base,
                        ch_types=None, n_chan=config.n_peak)

                    print('###\nPeak channels in Z-scored PSD for base'
                          'frequency %f: ' % basefreq)

                    # turn channel names into one list
                    # assume there was only one peak frequency
                    peak_ch_names_base = []
                    for chtype in peak_ch_types_base[0]:

                        peak_ch_names_base = peak_ch_names_base + peak_ch_types_base[0][chtype]

                    # Find channels with maximum Z-scores per channel type
                    # for oddball frequency
                    # "Latency" is frequency in Hz divided by 1000
                    peak_times_odd = [freq_odd]
                    peak_ch_types_odd = Ff.peak_channels_evoked(
                        evoked=evoked, peak_times=peak_times_odd,
                        ch_types=None, n_chan=config.n_peak)

                    print('###\nPeak channels in Z-scored PSD for oddball frequency %f: '
                          % freq_odd)

                    # turn channel names into one list
                    # assume there was only one peak frequency
                    peak_ch_names_odd = []
                    for chtype in peak_ch_types_odd[0]:

                        peak_ch_names_odd = peak_ch_names_odd + peak_ch_types_odd[0][chtype]

                    #

                    # Deepcopy because instance of evoked will be modified.
                    evoked = deepcopy(psd_z[fi])

                    # reduce evoked to peak channels for base frequency
                    evoked_peak = evoked.pick_channels(peak_ch_names_base)

                    psds[modal]['z_peak_base'][cond][freq].append(
                        evoked_peak)

                    #

                    # Deepcopy because instance of evoked will be modified.
                    evoked = deepcopy(psd_z[fi])

                    # reduce evoked to peak channels for oddball frequency
                    evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                    psds[modal]['z_peak_odd'][cond][freq].append(
                        evoked_peak)

                    #

                    evoked = deepcopy(psd[fi])

                    # base freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_base)

                    psds[modal]['peak_base'][cond][freq].append(
                        evoked_peak)

                    evoked = deepcopy(psd[fi])

                    # odd freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                    psds[modal]['peak_odd'][cond][freq].append(
                        evoked_peak)

                    #

                    evoked = deepcopy(psd_harm_odd[fi])

                    # base freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_base)

                    psds[modal]['harm_odd_peak_base'][cond][freq].append(
                        evoked_peak)

                    #

                    evoked = deepcopy(psd_harm_odd[fi])

                    # odd freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                    psds[modal]['harm_odd_peak_odd'][cond][freq].append(
                        evoked_peak)

                    #

                    evoked = deepcopy(psd_harm_base[fi])

                    # base freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_base)

                    psds[modal]['harm_base_peak_base'][cond][freq].append(
                        evoked_peak)

                    #

                    evoked = deepcopy(psd_harm_base[fi])

                    # odd freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                    psds[modal]['harm_base_peak_odd'][cond][freq].append(
                        evoked_peak)

                    #

                    evoked = deepcopy(psd_harm_topos_base[fi])

                    # base freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_base)

                    psds[modal]['peak_harm_topos_base'][cond][freq].append(
                        evoked_peak)

                    #

                    evoked = deepcopy(psd_harm_topos_odd[fi])

                    # base freq
                    evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                    psds[modal]['peak_harm_topos_odd'][cond][freq].append(
                        evoked_peak)

                    # Get data for group statistics (e.g. laterality)
                    # RMS across electrodes in ROI
                    for roi in config.electrode_ROIs:

                        ch_names = config.electrode_ROIs[roi]

                        for stim in ['base', 'odd']:

                            type_now = 'psd_harm_' + stim

                            evoked_roi = deepcopy(
                                psds[modal][type_now][cond][freq][-1])

                            evoked_roi.pick_channels(ch_names)

                            idx0 = evoked_roi.time_as_index(0.)

                            roi_rms = np.sqrt((evoked_roi.data[:, idx0]**2).mean())

                            roi_elecs[roi][cond][freq][stim][ss] = roi_rms

    # Reading source estimate (STC) data
    if 'stc' in modals:

        print('Reading source estimates.')

        modal = 'stc'

        for cond in conds:  # conditions

            print('Condition: %s.' % cond)

            if cond == 'face':  # no frequency sweep for faces

                freqs = ['6.0']  # base frequency for this condition (Hz, str)

                freq_odd = 1.2  # oddball frequency for this condition (Hz)

            else:  # for all word condition, use all sweep frequencies

                # base frequencies for this condition (Hz as string)
                freqs = freqs_all

                freq_odd = 1.0  # oddball frequency the same for all sweeps

            for (fi, freq) in enumerate(freqs):

                # for Evoked data are in one file for all frequencies
                # for STC data are in separate files per condition and freq
                for sbj_id in sbj_ids:  # across all subjects, EDIT ###

                    # path to subject's data
                    sbj_path = op.join(config.data_path,
                                       config.map_subjects[sbj_id][0])

                    print('Reading PSD results from STC files:')

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDTopo_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    psds[modal]['psd'][cond][freq].append(stc)

                    # fname_stc = op.join(
                    #     sbj_path, 'STC', '%sPSDTopoZ_%s_%s_mph-lh.stc' %
                    #     (prefix, cond, freq)
                    # )
                    # print(fname_stc)
                    # stc = mne.read_source_estimate(fname_stc)
                    # psds[modal]['psd_z'][cond][freq].append(stc)

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDHarm_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    psds[modal]['psd_harm_odd'][cond][freq].append(stc)

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDHarmBase_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    psds[modal]['psd_harm_base'][cond][freq].append(stc)

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumTopoOdd_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    psds[modal]['psd_sum_odd'][cond][freq].append(stc)

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumTopoBase_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    psds[modal]['psd_sum_base'][cond][freq].append(stc)

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumToposOdd_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    psds[modal]['psd_harm_topos_odd'][cond][freq].append(stc)

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumToposBase_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    psds[modal]['psd_harm_topos_base'][cond][freq].append(stc)

                # Grand-average STCs

                print('Grand-averaging source estimates.')

                for tt in stc_types:

                    stcs = psds[modal][tt][cond][freq]

                    avg_data = np.average([s.data for s in stcs], axis=0)

                    # turn average into source estimate object
                    stc_avg = SourceEstimate(avg_data, stcs[0].vertices,
                                             stcs[0].tmin, stcs[0].tstep)

                    fname_stc = op.join(config.grandmean_path, 'STC',
                                        '%s_%s_%s_%s' % (prefix, tt, cond,
                                                         freq))

                    print('Writing GM to %s.' % fname_stc)

                    stc_avg.save(fname_stc)

    # Compute Grand-Averages for Evoked data

    # Path for grand-mean results
    sbj_path = config.grandmean_path

    if 'evo_gm' in gm_modals:

        print('Grand-averaging evoked data.')

        psd_evo = psds['evo']

        for cond in conds:  # conditions

            print('###\nCondition: %s.\n###' % cond)

            if cond == 'face':  # no frequency sweep for faces

                freqs = ['6.0']  # base frequency for condition (Hz as str)

            else:  # for all word condition, use all sweep frequencies

                # base frequencies for this condition (Hz as string)
                freqs = freqs_all

            for tt in types:

                gm_evos = []  # get Evokeds for frequencies as list

                for freq in freqs:

                    # grand-average across subjects
                    evoked =\
                        mne.grand_average(psd_evo[tt][cond][freq])

                    evoked.comment = freq  # will be used in plotting script

                    # to keep everything
                    psds['evo_gm'][tt][cond][freq] = deepcopy(evoked)

                    # to write list of Evoked
                    gm_evos.append(evoked)

                fname_evo = op.join(sbj_path, 'AVE', 'GM_%s_%s-ave.fif' %
                                    (tt, cond))

                print('Writing GM to %s.' % fname_evo)

                mne.write_evokeds(fname=fname_evo, evoked=gm_evos)

            # the following cannot use grand_average() because channel
            # names differ across subjects
            # channel nambes can also differ across frequencies
            # therefore separate files for frequencies

            for tt in evo_types:

                for freq in freqs:

                    # Evokeds to average
                    evokeds = psd_evo[tt][cond][freq]

                    # grand-average across subjects
                    evoked =\
                        grand_average_evoked_arrays(evokeds)

                    fname_evo = op.join(
                        sbj_path, 'AVE', 'GM_%s_%s_%s-ave.fif' %
                        (tt, cond, freq))

                    print('Writing GM to %s.' % fname_evo)

                    mne.write_evokeds(fname=fname_evo, evoked=evoked)

            # Group Statistics for electrode ROIs
            for freq in freqs:

                for stim in ['base', 'odd']:

                    print('\nEEG laterality statistics for %s.' % stim)

                    data1 = roi_elecs['OT_L'][cond][freq][stim]
                    data2 = roi_elecs['OT_R'][cond][freq][stim]

                    stat, pv = ttest_rel(data1, data2)

                    print('T-test for L-R, %s | %s (%f vs %f).' %
                          (cond, freq, data1.mean(), data2.mean()))
                    print('%f, %f\n' % (stat, pv))

            # plot peak amplitudes across individual participants
            for freq in freqs:

                for stim in ['base', 'odd']:

                    type_now = 'harm_%s_peak_%s' % (stim, stim)

                    evokeds = deepcopy(psd_evo[type_now][cond][freq])

                    # get amplitudes at centre frequency per channel type
                    amps = get_amps_channel_types(evokeds)

                    for ch_type in amps.keys():

                        fig, ax = plt.subplots()

                        n = len(amps[ch_type])

                        x = np.arange(1, n + 1)

                        ax.bar(x, amps[ch_type])

                        threshold = 1.96
                        ax.plot([0., n], [threshold, threshold], "k--")

                        # output directory for figures
                        figs_path = op.join(config.grandmean_path, 'Figures_ICA')

                        fig_fname = op.join(
                            figs_path, 'face_amps_indiv_%s_%s.jpg' % (stim, ch_type))

                        print('Saving figure to %s.' % fig_fname)

                        fig.savefig(fig_fname)

                        plt.close(fig)

                    # put amplitudes into list of lists for correlation
                    amps_list = [amps['eeg'], amps['grad'], amps['mag']]

                    corrs = np.corrcoef(amps_list)

                    print('Condition: %s.' % stim)

                    print('Correlations of peak amplitudes between channel'
                          ' types across participants:')
                    print(corrs)



        # FOR FACES ONLY, put topographies for individual subjects together
        evos = psds['evo']['psd_sum_odd']['face']['6.0']

        data = evos[0].data

        # numpy array for topographies with shape (# sensors, # subjects)
        evo_mats = np.zeros((data.shape[0], len(evos)))

        for (ei, ee) in enumerate(evos):

            # evoked only has one sample
            evo_mats[:, ei] = ee.data[:, 0]

            evoked = EvokedArray(evo_mats, evos[0].info, tmin=0)

            fname_evo = op.join(
                sbj_path, 'AVE', 'GM_sum_indiv_topos_%s_%s-ave.fif' %
                ('face', '6.0'))

            print('Writing individual topographies to %s.' % fname_evo)

            mne.write_evokeds(fname_evo, evoked)

    return


def get_amps_channel_types(evokeds):
    """Extract RMS amplitudes across channels per channel type at latency 0.

    Parameters:
    evokeds: list of instances of Evoked
        Evokeds objects to extract amplitudes from, e.g. for peak channels
        per channel type.

    Returns:
    amps: dict of list
        Dictionary ('mag'/'grad'/'eeg') with lists (amplitudes per evoked).
    """

    ch_types = ['mag', 'grad', 'eeg']

    amps = {'mag': [], 'grad': [], 'eeg': []}

    for evoked in evokeds:

        for ch_type in ch_types:

            if ch_type in ['mag', 'grad']:

                meg = ch_type
                eeg = False

            else:

                eeg = True
                meg = False

            evo = deepcopy(evoked)

            evo.pick_types(meg=meg, eeg=eeg, eog=False, ecg=False)

            idx0 = evo.time_as_index(0.)

            rms = np.sqrt((evo.data[:, idx0]**2).mean())

            amps[ch_type].append(rms)

    return amps


def grand_average_conditions_data(data, evokeds, ch_names):
    """Average data arrays (e.g. peak channels) across subjects per condition.

    Parameters:
    data: dictionary of dictionary of lists of lists of numpy 2D arrays
        Dictionary contains conditions (e.g. "faces", "hflf"), then labels
        ('base'|'odd') then lists for all subjects and all sweep
        frequencies. The actual data are in numpy 2D (n_ch x n_t) arrays.
    evokeds: list of instances of Evoked
        Evokeds objects from which data was derived. For example, data may
        contain the data for peak channels extracted from evokeds.
    ch_names: list of str
        The instances of Evoked in evokes will be reduced to the channels
        specified in ch_names. The number of channels must be n_ch.

    Returns:
    gm_evokeds: dictionary of instances of Evoked
        Dictionary with conditions as in data, with data averaged across
        list items. The grand-averages are returned as
        gm_evoked[cond][sbj][freq].data where evoked is copied from last
        instance of evoked in evokeds per condition and frequency.
        Evokeds then contains appropriate info except for channel names (which
        will probably differ across subjects).
    """
    gm_evokeds = {}  # for instances of Evoked

    conds = data.keys()  # conditions in cond_evo

    # Slight complication: need to rearrange data in order to average

    datas = {}  # will collect data per condition, frequency, then subject
    for cond in conds:

        n_sbjs = len(data[cond])  # number of subjects

        # 'base'|'odd' depending on channel groups
        labels = data[cond].keys()

        gm_evokeds[cond] = {}  # Evokeds grand-averaged as Evoked
        datas[cond] = {}

        for lab in labels:

            gm_evokeds[cond][lab] = {}  # Evokeds grand-averaged as Evoked
            datas[cond][lab] = {}  # dict easier to handle below than list

        for ss in np.arange(n_sbjs):

            n_freqs = len(data[cond][list(labels)[0]][ss])

            for ff in np.arange(n_freqs):

                for lab in labels:

                    # data array for one frequency
                    data_freq = data[cond][lab][ss][ff]

                    if ss == 0:  # first time, initialise

                        datas[cond][lab][ff] = []

                    datas[cond][lab][ff].append(data_freq)

        # Grand-averaging across subjects
        for ff in np.arange(n_freqs):

            for lab in labels:

                # use last Evoked as template for this condition
                # the full evoked data with all channels doesn't contain labels
                # for channel groups
                evoked = deepcopy(evokeds[cond][-1][ff])

                # pick channel selection
                evoked.pick_channels(ch_names)

                # use as index, more informative later
                ff_str = evoked.comment

                # mean across list items
                gm_data = np.mean(datas[cond][lab][ff], axis=0)

                gm_evokeds[cond][lab][ff_str] = evoked

                # put peak channel data into Evoked object
                gm_evokeds[cond][lab][ff_str].data = gm_data

    return gm_evokeds


def grand_average_conditions_evo(evos):
    """Average evokeds across subjects per condition.

    Parameters:
    evos: dictionary of dictionary of lists
        Dictionary contains conditions (e.g. "faces", "hflf"), then lists that
        contain data for all subjects.

    Returns:
    gm_evo: dictionary
        Dictionary with conditions as in cond_evo, with data averaged across
        list items.
    """
    gm_evo = {}

    conds = evos.keys()  # conditions in cond_evo

    # Slight complication: need to rearrange Evoked object in order to use
    # mne.grand_average()

    evokeds = {}  # will collect Evoked per condition, frequency, then subject
    for cond in conds:

        n_sbjs = len(evos[cond])  # number of subjects

        evokeds[cond] = {}  # Evokeds per subject
        gm_evo[cond] = {}  # Evokeds grand-averaged

        for ss in np.arange(n_sbjs):

            n_freqs = len(evos[cond][ss])

            for ff in np.arange(n_freqs):

                # Evoked object for one frequency
                evo_freq = evos[cond][ss][ff]

                if ss == 0:  # first time, initialise
                    evokeds[cond][evo_freq.comment] = []

                evokeds[cond][evo_freq.comment].append(evo_freq)

        # Grand-averaging across subjects

        for ff in evokeds[cond].keys():

            gm_evo[cond][ff] = mne.grand_average(evokeds[cond][ff])

    return gm_evo


# get all input arguments except first
# if number not in config list, do all of them
if ((len(sys.argv) == 1) or
   (int(sys.argv[1]) > np.max(list(config.map_subjects.keys())))):

    # IDs don't start at 0
    sbj_ids = config.do_subjs

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


# requires all subjects to average across
grand_average_psds(sbj_ids)

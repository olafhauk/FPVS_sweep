
#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
Compute grand-average of the outputs of FPVS_PSD_sweep.py.

==========================================

OH, February 2020
"""

import sys

from os import path as op

import numpy as np

from matplotlib import pyplot as plt

from copy import deepcopy

from importlib import reload

import mne
from mne.report import Report

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


def grand_average_psds(sbj_ids):
    """Grand-average PSDs and derivatives across sbj_ids."""
    # initialise html report for one subject

    report = Report(subject='GM', title='FPVS PSDs GM')

    # get condition names and frequency names
    # assumed to be consistent across subjects
    sss_map_fname = config.sss_map_fnames[sbj_ids[0]]
    conds = []  # names of conditions
    for raw_stem_in in sss_map_fname[1][2:]:

        conds.append(raw_stem_in[:4])

    conds = np.unique(conds)

    # evoked files as dicts per condition, then lists across subjects
    psds_as_evo, psds_z_as_evo, sum_odd_as_evo, sum_base_as_evo,\
        psd_harm_as_evo, psd_harm_base_as_evo = {}, {}, {}, {}, {}, {}

    # data for peak channels per condition, then lists across subjects
    psds_data_peak, psds_z_data_peak, psd_harm_data_peak,\
        psd_harm_base_data_peak = {}, {}, {}, {}

    for cond in conds:  # conditions

        print('###\nCondition: %s.\n###' % cond)

        # create list of Evoked objects for all frequencies per condition
        # GM will be computed as average across list items
        psds_as_evo[cond], psds_z_as_evo[cond], sum_odd_as_evo[cond],\
            sum_base_as_evo[cond], psd_harm_as_evo[cond],\
            psd_harm_base_as_evo[cond] = [], [], [], [], [], []

        # The same but for arrays with peak data
        # for different channel groups
        psds_data_peak[cond] = {'base': [], 'odd': []}
        psds_z_data_peak[cond] = {'base': [], 'odd': []}
        psd_harm_data_peak[cond] = {'base': [], 'odd': []}
        psd_harm_base_data_peak[cond] = {'base': [], 'odd': []}

        if cond == 'face':  # no frequency sweep for faces

            freqs = ['6.0']  # base frequency for this condition (Hz as string)

            freq_odd = 1.2  # oddball frequency for this condition (Hz)

        else:  # for all word condition, use all sweep frequencies

            # base frequencies for this condition (Hz as string)
            freqs = freqs_all

            freq_odd = 1.0  # oddball frequency the same for all sweeps

        for sbj_id in sbj_ids:  # across all subjects, EDIT ###

            # path to subject's data
            sbj_path = op.join(config.data_path,
                               config.map_subjects[sbj_id][0])

            # Read Evoked objects for individual subjects:

            print('Reading PSD results from evoked files:')

            # separate filename prefixes for ICAed and non-ICAed data
            prefix = ''
            if 'ica' in config.raw_ICA_suff:
                prefix = 'ICA'

            fname_evo = op.join(sbj_path, '%sPSDTopo_%s%s' % (prefix,
                                                              cond, '-ave.fif'))
            print(fname_evo)
            psds_as_evo[cond].append(mne.read_evokeds(fname_evo))

            fname_evo = op.join(sbj_path, '%sPSDTopoZ_%s%s' % (prefix,
                                                               cond,
                                                               '-ave.fif'))
            print(fname_evo)
            psds_z_as_evo[cond].append(mne.read_evokeds(fname_evo))

            fname_evo = op.join(sbj_path, '%sPSDHarm_%s%s' % (prefix,
                                                              cond,
                                                              '-ave.fif'))
            print(fname_evo)
            psd_harm_as_evo[cond].append(mne.read_evokeds(fname_evo))

            fname_evo = op.join(sbj_path, '%sPSDHarmBase_%s%s' % (prefix,
                                                                  cond,
                                                                  '-ave.fif'))
            print(fname_evo)
            psd_harm_base_as_evo[cond].append(mne.read_evokeds(fname_evo))

            fname_evo = op.join(sbj_path, '%sPSDSumTopoOdd_%s%s' % (prefix,
                                                                    cond,
                                                                    '-ave.fif'))
            print(fname_evo)
            sum_odd_as_evo[cond].append(mne.read_evokeds(fname_evo))

            fname_evo = op.join(sbj_path, '%sPSDSumTopoBase_%s%s' % (prefix,
                                                                     cond,
                                                                     '-ave.fif'))
            print(fname_evo)
            sum_base_as_evo[cond].append(mne.read_evokeds(fname_evo))

            # check if there are as many evoked objects in list as their are
            # sweep frequencies for this condition
            if len(psds_as_evo[cond][-1]) != len(freqs):

                print('Number of evoked files (%d) and frequencies (%d) do '
                      'not match!' % (len(psds_as_evo), len(freqs)))
                return

            # get data for peak channels (may be different across subjects) for
            # later grand-averaging

            # data for peak channels per condition
            # dict for base and oddball frequeny, respectively
            # create list across base frequencies, then append per subject
            psds_tmp = {'base': [], 'odd': []}
            psds_z_tmp = {'base': [], 'odd': []}
            psd_harm_tmp = {'base': [], 'odd': []}
            psd_harm_base_tmp = {'base': [], 'odd': []}

            for (fi, freq) in enumerate(freqs):

                print('Frequency: %s.' % freq)

                # hack, float-to-string-to-float-again
                # to be consistent with FPVS_PSD_sweep_plot.py
                basefreq = float(freq)

                # Get max channels from z-scored PSD at base frequency
                # not oddball frequency, which would be biased.
                # This evoked is for condition cond, current subject -1 and
                # current frequency freq.
                evoked = deepcopy(psds_z_as_evo[cond][-1][fi])

                # Find channels with maximum Z-scores per channel type
                # for base frequency
                # "Latency" is frequency in Hz divided by 1000
                peak_times_base = [basefreq]
                peak_ch_types_base = Ff.peak_channels_evoked(evoked=evoked,
                                                             peak_times=peak_times_base,
                                                             ch_types=None,
                                                             n_chan=config.n_peak)

                print('###\nPeak channels in Z-scored PSD for base frequency %f: '
                      % basefreq)

                # turn channel names into one list
                # assume there was only one peak frequency
                peak_ch_names_base = []
                for chtype in peak_ch_types_base[0]:

                    peak_ch_names_base = peak_ch_names_base + peak_ch_types_base[0][chtype]

                # Find channels with maximum Z-scores per channel type
                # for oddball frequency
                # "Latency" is frequency in Hz divided by 1000
                peak_times_odd = [freq_odd]
                peak_ch_types_odd = Ff.peak_channels_evoked(evoked=evoked,
                                                            peak_times=peak_times_odd,
                                                            ch_types=None,
                                                            n_chan=config.n_peak)

                print('###\nPeak channels in Z-scored PSD for oddball frequency %f: '
                      % freq_odd)

                # turn channel names into one list
                # assume there was only one peak frequency
                peak_ch_names_odd = []
                for chtype in peak_ch_types_odd[0]:

                    peak_ch_names_odd = peak_ch_names_odd + peak_ch_types_odd[0][chtype]

                # Deepcopy because instance of evoked will be modified.
                evoked = deepcopy(psds_z_as_evo[cond][-1][fi])

                # reduce evoked to peak channels for base frequency
                evoked_peak = evoked.pick_channels(peak_ch_names_base)

                # get data as array for peak channels
                data_peak = evoked_peak.data

                psds_z_tmp['base'].append(data_peak)

                # Deepcopy because instance of evoked will be modified.
                evoked = deepcopy(psds_z_as_evo[cond][-1][fi])

                # reduce evoked to peak channels for oddball frequency
                evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                # get data as array for peak channels
                data_peak = evoked_peak.data

                psds_z_tmp['odd'].append(data_peak)

                #

                evoked = deepcopy(psds_as_evo[cond][-1][fi])

                # base freq
                evoked_peak = evoked.pick_channels(peak_ch_names_base)

                data_peak = evoked_peak.data

                psds_tmp['base'].append(data_peak)

                evoked = deepcopy(psds_as_evo[cond][-1][fi])

                # odd freq
                evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                data_peak = evoked_peak.data

                psds_tmp['odd'].append(data_peak)

                #

                evoked = deepcopy(psd_harm_as_evo[cond][-1][fi])

                # base freq
                evoked_peak = evoked.pick_channels(peak_ch_names_base)

                data_peak = evoked_peak.data

                psd_harm_tmp['base'].append(data_peak)

                evoked = deepcopy(psd_harm_as_evo[cond][-1][fi])

                # odd freq
                evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                data_peak = evoked_peak.data

                psd_harm_tmp['odd'].append(data_peak)

                #

                evoked = deepcopy(psd_harm_base_as_evo[cond][-1][fi])

                # base freq
                evoked_peak = evoked.pick_channels(peak_ch_names_base)

                data_peak = evoked_peak.data

                psd_harm_base_tmp['base'].append(data_peak)

                evoked = deepcopy(psd_harm_base_as_evo[cond][-1][fi])

                # odd freq
                evoked_peak = evoked.pick_channels(peak_ch_names_odd)

                data_peak = evoked_peak.data

                psd_harm_base_tmp['odd'].append(data_peak)

            for cc in psds_tmp.keys():  # for channel groups

                psds_data_peak[cond][cc].append(psds_tmp[cc])

                psds_z_data_peak[cond][cc].append(psds_z_tmp[cc])

                psd_harm_data_peak[cond][cc].append(psd_harm_tmp[cc])

                psd_harm_base_data_peak[cond][cc].append(psd_harm_base_tmp[cc])

    # GRAND-AVERAGE evoked data
    print('Computing Grand Averages for %d participants.' %
          len(psds_data_peak[cond]))

    psds_as_evo_gm = grand_average_conditions_evo(psds_as_evo)

    psds_z_as_evo_gm = grand_average_conditions_evo(psds_z_as_evo)

    sum_odd_as_evo_gm = grand_average_conditions_evo(sum_odd_as_evo)

    sum_base_as_evo_gm = grand_average_conditions_evo(sum_base_as_evo)

    psd_harm_as_evo_gm = grand_average_conditions_evo(psd_harm_as_evo)

    psd_harm_base_as_evo_gm =\
        grand_average_conditions_evo(psd_harm_base_as_evo)

    # GRAND-AVERAGE evoked data for peak channels
    psds_as_evo_peak_gm = grand_average_conditions_data(psds_data_peak,
                                                        psds_as_evo,
                                                        peak_ch_names_base)

    psds_z_as_evo_peak_gm = grand_average_conditions_data(psds_z_data_peak,
                                                          psds_z_as_evo,
                                                          peak_ch_names_base)

    psd_harm_as_evo_peak_gm = grand_average_conditions_data(psd_harm_data_peak,
                                                            psd_harm_as_evo,
                                                            peak_ch_names_base)

    psd_harm_base_as_evo_peak_gm =\
        grand_average_conditions_data(psd_harm_base_data_peak,
                                      psd_harm_base_as_evo,
                                      peak_ch_names_base)

    # Path for grand-mean results
    sbj_path = config.grandmean_path

    # output directory for figures
    figs_path = op.join(sbj_path, figs_dir)

    # PLOTTING ################################################################
    print('Plotting.')

    chtypes = ['mag', 'grad', 'eeg']

    for cond in conds:

        print('Condition %s.' % cond)

        # Plot topographies for sum across harmonic for oddball and base
        # frequencies
        for (odd_str, base_str) in zip(sum_odd_as_evo_gm[cond],
                                       sum_base_as_evo_gm[cond]):

            # topography for oddball frequency
            evoked = sum_odd_as_evo_gm[cond][odd_str]

            # get base frequency as string
            # remove possible '_' at beginning
            freq = odd_str[-4:].lstrip()  # '12.0', '_3.0' etc.

            if freq[0] == '_':

                freq = freq[1:]  # remove '_'

            print('freqstr, freq: %s, %s' % (odd_str, freq))

            # Use base frequency as "latency" in topoplot
            # Note: also for oddball frequency the "latency" is the base
            # frequency, because that's our experimental manipulation
            basefreq = float(freq)
            times = [basefreq]

            # remove '.'
            freq = ''.join(freq.split('.'))

            # label for condition and base frequency
            label_str = '%s_%s' % (cond, freq)

            file_label = odd_str

            # Filename stem for figure; channel type to be added later
            fname_fig = op.join(figs_path, prefix + file_label)

            print('Creating figure %s.' % fname_fig)

            figs = Ff.plot_evo_topomap(evoked, times, chtypes, fname_fig)

            for [fig, chtype] in zip(figs, chtypes):

                    sec_label = evoked.comment

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

            # topography for base frequency
            evoked = sum_base_as_evo_gm[cond][base_str]

            file_label = base_str

            # Filename stem for figure; channel type to be added later
            fname_fig = op.join(figs_path, prefix + file_label)

            print('Creating figure %s.' % fname_fig)

            figs = Ff.plot_evo_topomap(evoked, times, chtypes, fname_fig)

            for [fig, chtype] in zip(figs, chtypes):

                    sec_label = evoked.comment

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

        # plot evoked spectra and topographies (plot_joint())
        evo_list = [psds_as_evo_gm, psds_z_as_evo_gm]

        for evo in evo_list:

            for [fi, freqstr] in enumerate(evo[cond]):

                print('Plot GM evoked for %s %s.' % (cond, freqstr))

                evoked = evo[cond][freqstr]

                file_label = prefix + freqstr

                figs = Ff.plot_psd_as_evo(evoked, sbj_path, picks=None,
                                          txt_label=file_label,
                                          close_fig=close_fig,
                                          scalings=unit_scalings)

                for [fig, chtype] in zip(figs, chtypes):

                    sec_label = evoked.comment

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

                # plot spectra for EEG channel selections
                for roi in config.electrode_ROIs:

                    evoked_roi = deepcopy(evoked)

                    ch_names = config.electrode_ROIs[roi]

                    evoked_roi.pick_channels(ch_names)

                    # CROP PSD for display
                    evoked_roi.crop(tmin=config.crop_times[0],
                                    tmax=config.crop_times[1])

                    # Plot for peak channels without topographies
                    fig = evoked_roi.plot(spatial_colors=True, picks=None,
                                          scalings=unit_scalings, gfp=True)

                    fname_fig = op.join(figs_path, prefix +
                                        file_label + '_%s.pdf' % roi)

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    sec_label = evoked.comment + ' ' + roi

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

                plt.close('all')


        # plot evoked spectra for peak channels
        evo_list = [psds_as_evo_peak_gm, psds_z_as_evo_peak_gm]

        for evo in evo_list:

            # channel groups
            labels = evo[cond].keys()

            for lab in labels:

                for [fi, freqstr] in enumerate(evo[cond][lab]):

                    print('Plot GM evoked for %s %s %s.' % (cond, lab, freqstr))

                    evoked = evo[cond][lab][freqstr]

                    # CROP PSD for display
                    evoked.crop(tmin=config.crop_times[0],
                                tmax=config.crop_times[1])

                    # Plot for peak channels without topographies
                    fig = evoked.plot(spatial_colors=True, picks=None,
                                      scalings=unit_scalings, gfp=True)

                    filename = freqstr + '_Peak' + lab + '.pdf'

                    fname_fig = op.join(figs_path, filename)

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    sec_label = evoked.comment + ' Peak' + lab

                    report.add_figs_to_section(fig, sec_label, section=sec_label,
                                               scale=1)

            plt.close('all')

        # plot spectra around target frequencies
        evo_list = [psd_harm_as_evo_gm, psd_harm_base_as_evo_gm]

        for evo in evo_list:

            for [fi, freqstr] in enumerate(evo[cond]):

                print('Plot GM target frequencies for %s %s.' %
                      (cond, freqstr))

                evoked = evo[cond][freqstr]

                file_label = freqstr

                fig = evoked.plot(spatial_colors=True, picks=None,
                                  scalings=unit_scalings, gfp=True)

                fname_fig = op.join(figs_path, prefix + file_label + '.pdf')

                print('Creating figure %s.' % fname_fig)

                fig.savefig(fname_fig)

                sec_label = evoked.comment

                report.add_figs_to_section(fig, sec_label, section=sec_label,
                                           scale=1)

                # plot spectra for EEG channel selections
                for roi in config.electrode_ROIs:

                    evoked_roi = deepcopy(evoked)

                    ch_names = config.electrode_ROIs[roi]

                    evoked_roi.pick_channels(ch_names)

                    # Plot for peak channels without topographies
                    fig = evoked_roi.plot(spatial_colors=True, picks=None,
                                          scalings=unit_scalings, gfp=True)

                    fname_fig = op.join(figs_path, prefix + file_label +
                                        '_%s.pdf' % roi)

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    sec_label = evoked.comment + ' ' + roi

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

            plt.close('all')

        # plot spectra around target frequencies for peak channels
        evo_list = [psd_harm_as_evo_peak_gm, psd_harm_base_as_evo_peak_gm]

        for evo in evo_list:

            # channel groups
            labels = evo[cond].keys()

            for lab in labels:

                for [fi, freqstr] in enumerate(evo[cond][lab]):

                    print('Plot GM evoked for %s %s %s.' % (cond, lab, freqstr))

                    evoked = evo[cond][lab][freqstr]

                    # Plotting PSDs across harmonics
                    fig = evoked.plot(spatial_colors=True, picks=None,
                                      scalings=unit_scalings, gfp=True)

                    filename = prefix + freqstr + '_Peak' + lab + '.pdf'

                    fname_fig = op.join(figs_path, filename)

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    sec_label = evoked.comment + ' Peak' + lab
                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

                plt.close('all')

    # Save HTML report
    fname_report = op.join(figs_path, prefix + 'GM_report.html')

    report.save(fname_report, overwrite=True, open_browser=False)

    plt.close('all')

    return psds_as_evo_gm, psds_z_as_evo_gm, sum_odd_as_evo_gm,\
        sum_base_as_evo_gm, psd_harm_as_evo_gm, psd_harm_base_as_evo_gm


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
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


# requires all subjects to average across
psds_as_evo, psds_z_as_evo, sum_odd_as_evo, sum_base_as_evo, psd_harm_as_evo,\
    psd_harm_base_as_evo = grand_average_psds(sbj_ids)

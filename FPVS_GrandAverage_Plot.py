#!/imaging/local/software/miniconda/envs/mne0.21/bin/python
"""
Plot FPVS Grand-Mean data.
==========================================

OH, April 2020
"""

import os
from os import path as op

import numpy as np

import matplotlib
matplotlib.use('Agg') #  for running graphics on cluster
from matplotlib import pyplot as plt

# needed to run on SLURM
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from mayavi import mlab
mlab.options.offscreen = True

from copy import deepcopy

from importlib import reload

import mne
from mne.report import Report
from mne.source_estimate import SourceEstimate

import config_sweep as config
reload(config)

import FPVS_functions as Ff
reload(Ff)

print('Sunshine')

print(mne.__version__)

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

print(*freqs_all)

# separate filename prefixes for ICAed and non-ICAed data
prefix = ''
if 'ica' in config.raw_ICA_suff:
    prefix = 'ICA'

subjects_dir = config.subjects_dir

# average these three frequencies and plot separately
avg_freqs = ['6.0', '4.0', '3.0']

# output directory for figures
figs_path = op.join(config.grandmean_path, figs_dir)

# conditions
# conds = ['face', 'pwhf', 'pwlf', 'lfhf']
conds = config.do_conds

def grand_average_plot():
    """Plot grand-average PSDs and derivatives."""
    # initialise html report for one subject
    report = Report(subject='GM', title='FPVS PSDs GM')

    # for STC plotting
    subject = 'fsaverage'

    # # get condition names and frequency names from first subject
    # # assumed to be consistent across subjects
    # sss_map_fname = config.sss_map_fnames[1]
    # conds = []  # names of conditions
    # for raw_stem_in in sss_map_fname[1][2:]:

    #     conds.append(raw_stem_in[:4])

    # conds = np.unique(conds)

    # initialise

    # all psd results for evoked and STC
    # individual subjects and GM
    modals = ['evo', 'stc']
    gm_modals = ['evo_gm', 'stc_gm']
    # modals = ['stc']
    # gm_modals = ['stc_gm']

    # types = ['psd', 'psd_z', 'psd_sum_odd', 'psd_sum_base', 'psd_harm_odd',
    #          'psd_harm_base', 'psd_harm_topos_odd', 'psd_harm_topos_base']

    # evo_types = [
    #     'peak_odd', 'z_peak_odd', 'harm_odd_peak_odd', 'harm_base_peak_odd',
    #     'peak_base', 'z_peak_base', 'harm_odd_peak_base',
    #     'harm_base_peak_base', 'peak_harm_topos_odd', 'peak_harm_topos_base']

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

    psds = {}

    do_modals = modals + gm_modals

    # extract label amplitudes
    label_amps = {}
    for ss in stc_types:
        label_amps[ss] = {'lh': [], 'rh': []}

    # Initialise
    for modal in do_modals:

        psds[modal] = {}  # type of data

        do_types = types
        if modal[:3] == 'evo':  # add other types

            do_types = do_types + evo_types

        for tt in do_types:

            psds[modal][tt] = {}  # type of processed PSD

            for cond in conds:

                psds[modal][tt][cond] = {}  # sweep frequencies

                if cond == 'face':  # no frequency sweep for faces

                    # base frequency for this condition (Hz as string)
                    freqs = ['6.0']

                else:  # for all word condition, use all sweep frequencies

                    # base frequencies for this condition (Hz as string)
                    # add average across certain frequencies
                    freqs = freqs_all

                for freq in freqs:

                    psds[modal][tt][cond][freq] = []  # subjects

    # Read Evoked GM data

    # Path for grand-mean results
    sbj_path = config.grandmean_path

    if 'evo' in modals:

        modal = 'evo'  # do the evoked results here

        do_types = types + evo_types

        for tt in do_types:

            for cond in conds:  # conditions

                print('###\nCondition: %s.\n###' % cond)

                if cond == 'face':  # no frequency sweep for faces

                    # base frequency for this condition (Hz as string)
                    freqs = ['6.0']

                else:  # for all word condition, use all sweep frequencies

                    # base frequencies for this condition (Hz as string)
                    freqs = deepcopy(freqs_all)

                # if all frequencies in one evoked file
                if tt in types:

                    fname_evo = op.join(sbj_path, 'AVE', 'GM_%s_%s-ave.fif' %
                                        (tt, cond))

                    evokeds = mne.read_evokeds(fname=fname_evo)

                    for [fi, freq] in enumerate(freqs):

                        print(freq)

                        evoked = evokeds[fi]

                        print(evoked)

                        psds[modal][tt][cond][freq] = evoked

                elif tt in evo_types:

                    for [fi, freq] in enumerate(freqs):

                        fname_evo = op.join(
                            sbj_path, 'AVE', 'GM_%s_%s_%s-ave.fif' %
                            (tt, cond, freq))

                        evokeds = mne.read_evokeds(fname=fname_evo)

                        psds[modal][tt][cond][freq] = evokeds[0]

                print('Done reading evoked file.')

                # average certain frequencies, to be plotted separately
                if cond != 'face':  # if a word condition

                    print('Averaging frequencies: ')
                    print(*avg_freqs)

                    evo_freqs = []  # collect evoked across frequencies

                    for freq in avg_freqs:

                        # topography for oddball frequency
                        evoked = psds[modal][tt][cond][freq]

                        evo_freqs.append(evoked)

                    psds[modal][tt][cond]['avg'] =\
                        Ff.grand_average_evoked_arrays(evo_freqs)

        # PLOTTING ############################################################
        print('Plotting.')

        chtypes = ['mag', 'grad', 'eeg']  # for some topo plots

        # FOR FACES ONLY, put topographies for individual subjects together
        fname_evo = op.join(
            sbj_path, 'AVE', 'GM_sum_indiv_topos_%s_%s-ave.fif' %
            ('face', '6.0'))

        print('Reading evoked with topographies from %s.' % fname_evo)

        evoked = mne.read_evokeds(fname_evo, 0)

        print('Scaling topographies per sample.')
        evoked = Ff.scale_evoked_per_channel_type(evoked)

        for chtype in chtypes:

            # scaling to individual maxima per topography
            vmin, vmax = 0., 1.

            fig = evoked.plot_topomap(times=evoked.times, ch_type=chtype,
                                      vmin=vmin, vmax=vmax,
                                      scalings=unit_scalings[chtype],
                                      units='Z', show=False)

            fig_fname = op.join(
                figs_path, 'GM_sum_indiv_topos_%s_%s_%s.jpg' %
                ('face', '6.0', chtype))

            print('Saving individual topographies to %s.' % fig_fname)

            fig.savefig(fig_fname)

        # Plotting everything else

        for cond in conds:

            print('Condition %s.' % cond)

            # Plot topographies for sum across harmonic for oddball and base
            # frequencies

            do_types = ['psd_sum_odd', 'psd_sum_base']

            for tt in do_types:

                for freq in psds[modal][tt][cond]:

                    # topography
                    evoked = psds[modal][tt][cond][freq]

                    print('freq: %s' % str(freq))

                    times = [0.]

                    # remove '.'
                    freq_str = str(freq)
                    freq_str = ''.join(freq_str.split('.'))

                    sec_label = '%s_%s' % (cond, freq_str)

                    file_label = '%s_%s_%s_%s' % (prefix, cond, tt, freq_str)

                    # Filename stem for figure; channel type to be added later
                    fname_fig = op.join(figs_path, file_label)

                    print('Creating figure %s.' % fname_fig)

                    figs = Ff.plot_evo_topomap(evoked, times, chtypes,
                                               fname_fig)

                    for [fig, chtype] in zip(figs, chtypes):

                        report.add_figs_to_section(fig, tt, section=sec_label,
                                                   scale=1)

            # plot amplitudes across harmonics for electrode groups

            print('Plotting topographies and amplitudes across harmonics.')

            do_types = ['psd_harm_topos_base', 'psd_harm_topos_odd']

            for tt in do_types:

                for [fi, freq] in enumerate(psds[modal][tt][cond]):

                    # remove '.'
                    freq_str = str(freq)
                    freq_str = ''.join(freq_str.split('.'))

                    print('Plot GM evoked for %s %s.' % (cond, freq_str))

                    evoked = psds[modal][tt][cond][freq]

                    # change times for plotting to one sample per "second"
                    times = evoked.times
                    evoked.times = np.arange(0., len(times), 1.)

                    # label for condition and base frequency
                    sec_label = '%s_%s' % (cond, freq_str)

                    file_label = '%s_%s_%s_%s' % (prefix, cond, tt, freq_str)

                    # Plot topopraphies for all harmonics

                    # Filename stem for figure; channel type to be added later
                    fname_fig = op.join(figs_path, file_label)

                    print('Creating figure %s.' % fname_fig)

                    times = evoked.times  # all harmonics

                    figs = Ff.plot_evo_topomap(evoked, times, chtypes,
                                               fname_fig)

                    # plot spectra for EEG channel selections
                    for roi in config.electrode_ROIs:

                        evoked_roi = deepcopy(evoked)

                        ch_names = config.electrode_ROIs[roi]

                        evoked_roi.pick_channels(ch_names)

                        # Plot for peak channels without topographies
                        fig = evoked_roi.plot(spatial_colors=True, picks=None,
                                              scalings=unit_scalings,
                                              gfp=False)

                        fname_fig = op.join(figs_path,
                                            file_label + '_%s.jpg' % roi)

                        print('Creating figure %s.' % fname_fig)

                        fig.savefig(fname_fig)

                        sec_label = sec_label + ' ' + roi

                        report.add_figs_to_section(fig, sec_label,
                                                   section=sec_label, scale=1)

                    # get singular values per channel type
                    # don't include last harmonic because of MEG artefact
                    idx = np.arange(0, evoked.data.shape[1] - 1, 1)
                    ss = Ff.svd_per_channel_type(evoked, idx)[0]

                    # channel types for SVD
                    ch_types = ['grad', 'mag', 'eeg']

                    # create new pyplot figure, subplots for channel types
                    fig, axs = plt.subplots(len(ch_types), 1)

                    for [ci, cht] in enumerate(ch_types):

                        # turn singular values into variances
                        s = 100. * ss[cht]**2 / (ss[cht]**2).sum()

                        x = np.arange(1, len(s) + 1, 1)

                        # plot singular values to figure
                        axs[ci].plot(x, s)

                        axs[ci].set_title(cht)

                    fig.tight_layout(pad=1.)

                    fname_fig = op.join(
                        figs_path, file_label + '_svd.jpg')

                    # save figure for this channel type
                    fig.savefig(fname_fig)

                    plt.close('all')  # close pyplot figures

            # plot evoked spectra and topographies (plot_joint())
            do_types = ['psd', 'psd_z']

            for tt in do_types:

                for [fi, freq] in enumerate(psds[modal][tt][cond]):

                    # remove '.'
                    freq_str = str(freq)
                    freq_str = ''.join(freq_str.split('.'))

                    print('Plot GM evoked for %s %s.' % (cond, freq_str))

                    evoked = psds[modal][tt][cond][freq]

                    # label for condition and base frequency
                    sec_label = '%s_%s' % (cond, freq_str)

                    file_label = '%s_%s_%s_%s' % (prefix, cond, tt, freq_str)

                    figs = Ff.plot_psd_as_evo(evoked, sbj_path, picks=None,
                                              txt_label=file_label,
                                              close_fig=close_fig,
                                              scalings=unit_scalings)

                    for [fig, chtype] in zip(figs, chtypes):

                        report.add_figs_to_section(fig, file_label,
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
                                              scalings=unit_scalings,
                                              gfp=False)

                        fname_fig = op.join(figs_path,
                                            file_label + '_%s.jpg' % roi)

                        print('Creating figure %s.' % fname_fig)

                        fig.savefig(fname_fig)

                        sec_label = sec_label + ' ' + roi

                        report.add_figs_to_section(fig, sec_label,
                                                   section=sec_label, scale=1)

                    plt.close('all')

            # plot evoked spectra for peak channels
            do_types = ['peak_odd', 'peak_base', 'z_peak_odd', 'z_peak_base']

            for tt in do_types:

                for [fi, freq] in enumerate(psds[modal][tt][cond]):

                    # remove '.'
                    freq_str = str(freq)
                    freq_str = ''.join(freq_str.split('.'))

                    print('Plot GM evoked for %s %s.' % (cond, freq_str))

                    evoked = psds[modal][tt][cond][freq]

                    # CROP PSD for display
                    evoked.crop(tmin=config.crop_times[0],
                                tmax=config.crop_times[1])

                    # Plot for peak channels without topographies
                    fig = evoked.plot(spatial_colors=True, picks=None,
                                      scalings=unit_scalings, gfp=False)

                    sec_label = '%s_%s' % (cond, freq_str)

                    file_label = '%s_%s_%s_%s' % (prefix, cond, tt, freq_str)

                    fname_fig = op.join(figs_path, file_label + '.jpg')

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

            plt.close('all')

            # plot amplitudes of harmonics for peak channels
            do_types = ['peak_harm_topos_odd', 'peak_harm_topos_base']

            for tt in do_types:

                for [fi, freq] in enumerate(psds[modal][tt][cond]):

                    # remove '.'
                    freq_str = str(freq)
                    freq_str = ''.join(freq_str.split('.'))

                    print('Plot GM evoked for %s %s.' % (cond, freq_str))

                    evoked = psds[modal][tt][cond][freq]

                    times = evoked.times
                    evoked.times = np.arange(0., len(times), 1.)

                    # Plot for peak channels without topographies
                    fig = evoked.plot(spatial_colors=True, picks=None,
                                      scalings=unit_scalings, gfp=False,
                                      sphere=0.)

                    sec_label = '%s_%s' % (cond, freq_str)

                    file_label = '%s_%s_%s_%s' % (prefix, cond, tt, freq_str)

                    fname_fig = op.join(figs_path, file_label + '.jpg')

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    # also create PDF because some edits may be needed
                    fname_fig = op.join(figs_path, file_label + '.pdf')

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

            plt.close('all')

            # plot spectra around target frequencies
            do_types = ['psd_harm_odd', 'psd_harm_base']

            for tt in do_types:

                for [fi, freq] in enumerate(psds[modal][tt][cond]):

                    # remove '.'
                    freq_str = str(freq)
                    freq_str = ''.join(freq_str.split('.'))

                    print('Plot GM target frequencies for %s %s.' %
                          (cond, freq_str))

                    evoked = psds[modal][tt][cond][freq]

                    fig = evoked.plot(spatial_colors=True, picks=None,
                                      scalings=unit_scalings, gfp=False)

                    sec_label = '%s_%s' % (cond, freq_str)

                    file_label = '%s_%s_%s_%s' % (prefix, cond, tt, freq_str)

                    fname_fig = op.join(figs_path, file_label + '.jpg')

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

                    # plot spectra for EEG channel selections
                    for roi in config.electrode_ROIs:

                        evoked_roi = deepcopy(evoked)

                        ch_names = config.electrode_ROIs[roi]

                        evoked_roi.pick_channels(ch_names)

                        # Plot for peak channels without topographies
                        fig = evoked_roi.plot(spatial_colors=True, picks=None,
                                              scalings=unit_scalings,
                                              gfp=False)

                        fname_fig = op.join(figs_path,
                                            file_label + '_%s.jpg' % roi)

                        print('Creating figure %s.' % fname_fig)

                        fig.savefig(fname_fig)

                        sec_label = sec_label + ' ' + roi

                        report.add_figs_to_section(fig, sec_label,
                                                   section=sec_label, scale=1)

                plt.close('all')

            # plot spectra around target frequencies for peak channels
            do_types = ['harm_odd_peak_odd', 'harm_base_peak_odd',
                        'harm_odd_peak_base', 'harm_base_peak_base']

            for tt in do_types:

                for [fi, freq] in enumerate(psds[modal][tt][cond]):

                    # remove '.'
                    freq_str = str(freq)
                    freq_str = ''.join(freq_str.split('.'))

                    print('Plot GM evoked for %s %s.' % (cond, freq_str))

                    evoked = psds[modal][tt][cond][freq]

                    # Plotting PSDs across harmonics
                    fig = evoked.plot(spatial_colors=True, picks=None,
                                      scalings=unit_scalings, gfp=False)

                    sec_label = '%s_%s' % (cond, freq_str)

                    file_label = '%s_%s_%s_%s' % (prefix, cond, tt,
                                                  freq_str)

                    fname_fig = op.join(figs_path, file_label + '.jpg')

                    print('Creating figure %s.' % fname_fig)

                    fig.savefig(fname_fig)

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

                plt.close('all')

    # Plot STCs

    if 'stc' in modals:

        modal = 'stc'  # do source estimates here

        for tt in stc_types:

            for cond in conds:  # conditions

                print('###\nCondition: %s.\n###' % cond)

                if cond == 'face':  # no frequency sweep for faces

                    freqs = ['6.0']  # base frequency for condition (Hz as str)

                else:  # for all word condition, use all sweep frequencies

                    # base frequencies for this condition (Hz as string)
                    freqs = deepcopy(freqs_all)

                stc_freqs = {}  # collect STCs for all frequencies

                for (fi, freq) in enumerate(freqs):

                    fname_stc = op.join(
                        config.grandmean_path, 'STC',
                        '%s_%s_%s_%s-lh.stc' % (prefix, tt, cond, freq)
                    )

                    print('Reading source estimate from %s.' % fname_stc)

                    stc = mne.read_source_estimate(fname_stc)

                    stc_freqs[freq] = stc

                # for word conditions average certain frequencies
                if cond != 'face':

                    # pick STCs for certain frequencies for averaging
                    stcs = [stc_freqs[ff] for ff in avg_freqs]

                    avg_data = np.average([s.data for s in stcs], axis=0)

                    # turn average into source estimate object
                    stc_freqs['avg'] = SourceEstimate(
                        avg_data, stcs[0].vertices, stcs[0].tmin,
                        stcs[0].tstep)

                    # include average frequency from now on
                    freqs.append('avg')

                for (fi, freq) in enumerate(freqs):

                    # use STC for this frequency
                    stc = stc_freqs[freq]

                    time_label = None  # '%s %s' % (cond, freq)

                    # index to time point 0, which will be plotted
                    idx0 = np.abs(stc.times).argmin()

                    thresh = np.abs(stc.data[:, idx0]).max()

                    # # get some round numbers for colour bar

                    # if thresh < 10:

                    #     thresh = np.floor(thresh)

                    # elif thresh < 50:

                    #     thresh = 5 * np.floor(thresh / 5.)

                    # else:

                    #     thresh = 10 * np.floor(thresh / 10.)

                    # plot for left and right hemisphere
                    for hemi in ['both']:  # ['lh', 'rh']:

                        # for some reason, 'both' only works for 'ven' but not
                        # for 'lat'
                        for view in ['ven']:

                            brain = stc.plot(
                                subject=subject, initial_time=0.,
                                time_label=time_label,
                                subjects_dir=subjects_dir,
                                clim=dict(kind='value',
                                          lims=[0, thresh / 2., thresh]),
                                hemi=hemi, views=view
                            )

                            fname_fig = op.join(
                                figs_path,
                                '%s_%s_%s_%s_STC_%s_%s.jpg' %
                                (prefix, tt, cond, freq, hemi, view)
                            )

                            print('Saving figure to %s.' % fname_fig)

                            mlab.savefig(fname_fig)

                            mlab.close(all=True)

                    # plot for left and right hemisphere
                    for hemi in ['lh', 'rh']:

                        # for some reason, 'both' only works for 'ven' but not
                        # for 'lat'
                        for view in ['lat']:

                            # apparently 'brain' required for saving?
                            brain = stc.plot(
                                subject=subject, initial_time=0.,
                                time_label=time_label,
                                subjects_dir=subjects_dir,
                                clim=dict(kind='value',
                                          lims=[0, thresh / 2., thresh]),
                                hemi=hemi, views=view
                            )

                            fname_fig = op.join(
                                figs_path,
                                '%s_%s_%s_%s_STC_%s_%s.jpg' %
                                (prefix, tt, cond, freq, hemi, view)
                            )

                            print('Saving figure to %s.' % fname_fig)

                            mlab.savefig(fname_fig)

                            mlab.close(all=True)

    # Save HTML report
    fname_report = op.join(figs_path, prefix + 'GM_report.html')

    report.save(fname_report, overwrite=True, open_browser=False)

    plt.close('all')

    return

grand_average_plot()

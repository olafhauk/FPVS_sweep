#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
Compute PSD for average raw data for FPVS Frequency Sweep.

Average raw data from FPVS_get_sweeps.py.
Plot figures.
Compute z-scores.
Compute TFR if specified.
==========================================

OH, October 2019
"""

import sys

import os
from os import path as op
import numpy as np

from copy import deepcopy

# import matplotlib
# matplotlib.use('Agg') #  for running graphics on cluster ### EDIT

from matplotlib import pyplot as plt

from importlib import reload

import mne
from mne.report import Report

import config_sweep as config
reload(config)


print(mne.__version__)

# perform TFR of raw data or not
# do_tfr = config.do_tfr

# sub-directory for figures per subject
figs_dir = 'Figures_ICA'
# figs_dir = 'Figures'

close_fig = 1  # close figures only if close_fig==1

# plt.ion() # interactive plotting

# for some plots of SNRs
unit_scalings = dict(eeg=1., mag=1., grad=1.)

def run_PSD_raw(sbj_id):
    """Compute spectra for one subject."""

    # initialise html report for one subject
    report = Report(subject=str(sbj_id), title='FPVS PSDs')

    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    # get condition names and frequency names
    conds = []  # names of conditions
    for raw_stem_in in sss_map_fname[1][2:]:

        conds.append(raw_stem_in[:4])

    conds = np.unique(conds)

    freqs_all = [str(ff) for ff in config.fpvs_freqs]

    print(freqs_all)

    # initialise sum across harmonics for conditions
    sum_harms = {}
    for cond in conds:

        sum_harms[cond] = {}

    # Go through conditions and frequencies
    for cond in conds:  # conditions

        if cond == 'face':  # hack, no frequency sweep for faces

            freqs = ['6.0']

        else:  # for all word condition, use all sweep frequencies

            freqs = freqs_all

        for freq in freqs:  # frequencies

            sum_harms[cond][freq] = []  # initialise for this base frequency

            # remove dot from frequency string

            fname = 'rawavg_%s_%s.fif' % (cond, ''.join(freq.split('.')))
            fname_raw_in = op.join(sbj_path, fname)

            print('Reading average raw data from %s:' % fname_raw_in)

            raw = mne.io.read_raw_fif(fname_raw_in, preload=True)

            # picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False,
            #                       ecg=False, stim=False)

            # reduce raw data to relevant channels
            raw.pick_types(meg=True, eeg=True, eog=False, ecg=False,
                           stim=False, misc=False, chpi=False)

            # Compute PSD for raw data
            n_fft = len(raw.times)

            print('Computing psd_welch().')
            psds, psd_freqs = mne.time_frequency.psd_welch(raw, fmin=0.,
                                                           fmax=40.,
                                                           n_fft=n_fft)

            print('Frequency resolution: %f.' % (psd_freqs[1] - psd_freqs[0]))

            info = raw.info
            # pretend sample frequency for frequencies
            info['sfreq'] = 1000. / (psd_freqs[1] - psd_freqs[0])

            # convert to amplitudes (rather than power)
            psds = np.sqrt(psds)

            # "BASELINE-CORRECT" PSDs with neighbouring frequency bins
            print("\nCorrecting baseline.\n")
            if config.psd_base_bins > 0:

                psds = psd_correct_baseline(psds, config.psd_base_bins,
                                            config.psd_n_gap)

            print('Computing SNRs.')
            psds_z = psd_convert_to_snr(psds, config.psd_snr_bins,
                                        config.psd_n_gap)

            # Raw PSDs
            psds_as_evo = mne.EvokedArray(psds, info,
                                          tmin=psd_freqs[0] / 1000.,
                                          comment=('PSD ' + cond + ' ' + freq))
            # z-scored PSDs
            psds_z_as_evo = mne.EvokedArray(psds_z, info,
                                            tmin=psd_freqs[0] / 1000.,
                                            comment=('PSD_Z ' + cond + ' ' +
                                                     freq))

            # Plot PSD as spectrum plus topographies (plot_joint())
            print('Plotting PSDs.')

            # label for condition and base frequency
            label_str = '%s_%s' % (cond, ''.join(freq.split('.')))

            file_label = 'PSDtopo_%s' % label_str
            figs = plot_psd_as_evo(psds_as_evo, sbj_path, file_label,
                                   close_fig)

            chtypes = ['mag', 'grad']  # sequence of plots in plot_joint
            if len(figs) == 3:
                chtypes = chtypes + ['eeg']

            for [fig, chtype] in zip(figs, chtypes):

                sec_label = '%s_%s' % (label_str, chtype)

                report.add_figs_to_section(fig, sec_label, section=sec_label,
                                           scale=1)

            file_label = 'PSDtopoZ_%s' % label_str
            figs = plot_psd_as_evo(psds_z_as_evo, sbj_path, file_label,
                                   close_fig, scalings=unit_scalings)

            for [fig, chtype] in zip(figs, chtypes):

                sec_label = '%s_%s_Z' % (label_str, chtype)

                report.add_figs_to_section(fig, sec_label, section=sec_label,
                                           scale=1)

            # Compute the sum across harmonics of oddball frequency for this
            # condition and base frequency

            print('Combining harmonics.')

            if cond == 'face':

                # round to make sure combine_harmonics finds the right
                # frequencies
                oddfreq = round(config.fpvs_odd_freq['faces'], 2)

            else:

                oddfreq = round(config.fpvs_odd_freq['words'], 2)

            basefreq = float(freq)  # hack, float-to-string-to-float-again

            sum_harms[cond][freq] = combine_harmonics(psds=psds_z,
                                                      freqs=psd_freqs,
                                                      basefreq=basefreq,
                                                      oddfreq=oddfreq,
                                                      n_harms=config.fpvs_n_harms,
                                                      method='sum')

            to_evo = np.expand_dims(sum_harms[cond][freq], 1)

            sum_as_evo = mne.EvokedArray(to_evo, info,
                                         tmin=basefreq / 1000.,
                                         comment=('Sum ' + cond + ' ' + freq))

            # get PSDs around harmonics
            psd_harm = psds_for_harmonics(psds=psds_z, freqs=psd_freqs,
                                          basefreq=basefreq, oddfreq=oddfreq,
                                          n_harms=config.fpvs_n_harms,
                                          n_bins=config.psd_snr_bins,
                                          method='sum')

            # for testing, do it for base frequency
            psd_harm_base = psds_for_harmonics(psds=psds_z, freqs=psd_freqs,
                                               basefreq=999., oddfreq=basefreq,
                                               n_harms=3,
                                               n_bins=config.psd_snr_bins,
                                               method='sum')

            tmin = -config.psd_snr_bins * 0.001
            info['sfreq'] = 1000.  # to display samples as time points
            psd_harm_as_evo = mne.EvokedArray(psd_harm, info, tmin=tmin,
                                              comment=('PSD Harm ' + cond + ' ' + freq))

            psd_harm_base_as_evo = mne.EvokedArray(psd_harm_base, info, tmin=tmin,
                                                   comment=('PSD Harm ' + cond + ' ' + freq))

            # Plotting PSDs across harmonics
            fig = psd_harm_as_evo.plot(spatial_colors=True,
                                       scalings=unit_scalings)

            # path to sub-directory for figures
            figs_path = op.join(sbj_path, figs_dir)

            fname_fig = op.join(figs_path, 'HarmPSD_%s_%s.pdf'
                                % (cond, freq))

            fig.savefig(fname_fig)

            sec_label = '%s_%s_HarmPSD' % (label_str, chtype)
            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            plt.close(fig)

            # Plotting PSDs across harmonics
            fig = psd_harm_base_as_evo.plot(spatial_colors=True,
                                            scalings=unit_scalings)

            fname_fig = op.join(figs_path, 'HarmPSDbase_%s_%s.pdf'
                                % (cond, freq))

            fig.savefig(fname_fig)

            sec_label = '%s_%s_HarmPSDbase' % (label_str, chtype)
            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            plt.close(fig)

            times = [basefreq / 1000.]

            for chtype in chtypes:

                fig = sum_as_evo.plot_topomap(times=times, ch_type=chtype,
                                              scalings=unit_scalings[chtype],
                                              units='snr',
                                              show=False)

                if not os.path.isdir(figs_path):
                    os.mkdir(figs_path)

                fname_fig = op.join(figs_path, 'SumTopo_%s_%s_%s.pdf'
                                    % (cond, freq, chtype))

                print('Saving combined topography across harmonics to: %s.' %
                      fname_fig)

                fig.savefig(fname_fig)

                sec_label = '%s_%s_Sum' % (label_str, chtype)

                report.add_figs_to_section(fig, sec_label, section=sec_label,
                                           scale=1)

                plt.close(fig)

        # Save HTML report
        fname_report = op.join(figs_path, str(sbj_id) +
                               '_report.html')

        report.save(fname_report, overwrite=True, open_browser=False)

    return


    #     # find first trigger at beginning of FPVS sequence
    #     start_idx = np.where(events[:, 2] == 19)[0]

    #     # find trigger at end of FPVS sequence
    #     end_idx = np.where(events[:, 2] == 20)[0]

    #     # time between first trigger and first stimulus onset in seconds
    #     delay = (events[start_idx + 1, 0] - events[start_idx, 0]) / raw.info['sfreq']

    #     # ignore first cycle of stimulus sequence
    #     # different for pictures and words
    #     if raw_stem_in[0] == 'w':

    #         delay = delay + 2.

    #     else:

    #         delay = delay + 1.5

    #     print('Initial delay: %s' % delay)

    #     # start and end of sequence in seconds
    #     # 50s window, integer of presentation cycle (see Lochy et al., Retter&Rossion 2016)
    #     # takes into account initial delay
    #     start_lat = float(((events[start_idx, 0] - raw.first_samp) / raw.info['sfreq']) + delay )

    #     end_lat = start_lat + 50.

    #     # end_lat = ((events[end_idx,0] - raw.first_samp) / raw.info['sfreq']) - 2.

    #     # length of raw data segment in seconds
    #     duration = end_lat - start_lat

    #     # number of samples in raw data segment
    #     # n_fft = int(duration * raw.info['sfreq'])
    #     n_fft = 50000  # 2**15 # closest to 50000

    #     print('Cropping raw data between %f and %f (duration %f).' % (start_lat, end_lat, duration))

    #     raw.crop(start_lat, end_lat)

    #     # PSD subplots for channel types separately

    #     # Pick channel selection for MEG
    #     # note: for PSD important to exclude EOG etc.
    #     channel_groups = ['Left-occipital', 'Right-occipital', 'Left-temporal', 'Right-temporal']
    #     selection_meg = []
    #     for chg in channel_groups:

    #         selection_meg = selection_meg + mne.read_selection(chg)

    #     # remove blanks from channel names
    #     selection_meg = [ss.replace(' ', '') for ss in selection_meg]

    #     # number of subplots depends on presence of EEG
    #     ch_types = ['mag', 'grad']

    #     if is_eeg:
    #         ch_types = ch_types + ['EEG']

    #     # create subplot with one pannel per channel type
    #     # fig[raw_stem_in], ax = plt.subplots(len(ch_types),1,sharex=True)

    #     # # Fill subplots
    #     # for (chi,cht) in enumerate(ch_types):

    #     #     # choose flags for "picks"
    #     #     if cht=='EEG':
    #     #         eeg = True
    #     #         meg = False
    #     #         selection = config.occ_eeg
    #     #     else:
    #     #         eeg=False
    #     #         meg = cht
    #     #         selection = selection_meg

    #     #     picks = mne.pick_types(raw.info, meg=meg, eeg=eeg, eog=False, ecg=False,
    #     #                                 stim=False, exclude='bads', selection=selection)

    #     #     raw.plot_psd(fmin=0.5, fmax=11., average=False, dB=False, ax=ax[chi],
    #     #                                        n_fft=n_fft, picks=picks)

    #     # picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False, ecg=False, stim=False)

    #     # reduce raw data to relevant channels
    #     raw.pick_types(meg=meg, eeg=eeg, eog=False, ecg=False, stim=False, misc=False, chpi=False)

    #     # Compute PSD for raw data
    #     psds, freqs = mne.time_frequency.psd_welch(raw, fmin=0., fmax=11., n_fft=n_fft)

    #     # Turn PSD into Evoked object for more visualisation options (like topographies)
        
    #     # ch_names = raw.info['ch_names']
    #     # ch_types = raw.info['ch_types']
    #     # sfreq = 1000.*(freqs[1]-freqs-[0]) # pretend sample frequency for frequencies
    #     # info = mne.create_info(ch_names, sfreq, ch_types) # "ch_types" a problem

    #     info = raw.info
    #     info['sfreq'] = 1000./(freqs[1]-freqs[0]) # pretend sample frequency for frequencies

    #     # convert to amplitudes (rather than power)
    #     psds = np.sqrt(psds)

    #     # "BASELINE-CORRECT" PSDs with neighbouring frequency bins
    #     print("\nCorrecting baseline.\n")
    #     if config.psd_base_bins > 0:

    #         psds = psd_correct_baseline(psds)

    #     # convert PSD to Evoked object
    #     info_copy = deepcopy(info) # to change sampling frequency
    #     info_copy['sfreq'] = 1000./(freqs[1]-freqs[0])
    #     psds_as_evo = mne.EvokedArray(psds, info_copy, tmin=freqs[0]/1000., comment='PSD '+ raw_stem_in)

        
    #     ## built-in plots for raw PSDs don't plot topographies

    #     # PSDs per sensor (for some reason flat line for dB=False)
    #     # fig[raw_stem_in] = raw.plot_psd_topo(fmin=0., fmax=12., dB=True, n_fft=n_fft)

    #     # fig[raw_stem_in].suptitle(raw_stem_in)

    #     # # Plot PSDs
    #     # fig = psds_as_evo.plot(spatial_colors=True, units=units, window_title=raw_stem_in)

    #     # fname_fig = op.join(sbj_path, 'Figures', 'PSDs_' + raw_stem_in + '_sss_f_raw.pdf')
    #     # print('Saving figure to %s' % fname_fig)

    #     # fig.savefig(fname_fig)        

    #     # # Plot topographies for channel types separately
    #     # fig, ax = plt.subplots(len(units.keys()),len(times),sharex=True)

    #     # for (ci,ch_type) in enumerate(units.keys()):

    #     #     psds_as_evo.plot_topomap(times=times, ch_type=ch_type, scalings=scalings, colorbar=False,
    #     #                              time_format='%.2f Hz', axes=ax[ci], title=raw_stem_in)

    #     # Plot PSD as spectrum plus topographies (plot_joint())
    #     plot_psd_as_evo(psds_as_evo, sbj_path, raw_stem_in, '', close_fig)

    #     # # keep a copy for scaling below
    #     # psds_tmp = deepcopy(psds_as_evo)

    #     # # CROP PSD for display
    #     # psds_as_evo.crop(tmin=config.crop_times[0], tmax=config.crop_times[1])

    #     # # plot topographic maps for different frequencies
    #     # data = psds_as_evo.data

    #     # # quick hack to get scaling approximately right
    #     # # avoid first 1Hz in scaling
    #     # psds_tmp.crop(tmin=0.001, tmax=config.crop_times[1])

    #     # # Default mne-python scalings for evo.plot(), just for clarity
    #     # scalings = dict(eeg=1e6, grad=1e13, mag=1e15)
    #     # # show y-axis with original values
    #     # # scalings = dict(eeg=1., grad=1., mag=1.)

    #     # # plot y-axis range
    #     # ylim = {'mag': [0,np.max(psds_tmp.data[2:306:3,:])*scalings['mag']],
    #     #         'grad': [0,np.max([psds_tmp.data[0:306:3,:], psds_tmp.data[1:306:3,:]])*scalings['grad']]}
    #     # if is_eeg:
    #     #     ylim['eeg'] = [0,np.max(psds_tmp.data[306:376,:]*scalings['eeg'])]

    #     # print('Upper limits for Mag: %e, Grad: %e ' % (ylim['mag'][1], ylim['grad'][1]))
    #     # if is_eeg:
    #     #     print('EEG: %e.\n' % ylim['eeg'][1])

    #     # # ylim = {}
    #     # # if is_eeg:
    #     # #     ylim['eeg'] = [0, 1e-9]
        
    #     # # ylim['mag'] = [0, 1e-9]

    #     # # ylim['grad'] = [0, 1e-9]

    #     # ts_args = dict(spatial_colors=True, scalings=scalings, units=units, ylim=ylim)
    #     # topomap_args = dict(scalings=scalings, time_format='%.2f Hz')

    #     # fig = psds_as_evo.plot_joint(times=config.times, title=raw_stem_in, ts_args=ts_args, topomap_args=topomap_args)

    #     # crop_str = str(int(10000.*config.crop_times[0])) + '_' + str(int(10000.*config.crop_times[1]))
    #     # ch_types = ['mag', 'grad'] # apparent sequence of plots in plot_joint output
    #     # if len(fig)==3:
    #     #     ch_types = ['eeg'] + ch_types

    #     # for (cc,ff) in zip(ch_types, fig):

    #     #     for fig_format in ['.pdf']: # jpg doesn't work, png does

    #     #         # name depends on crop_time in config.py
    #     #         fname_fig = op.join(sbj_path, 'Figures', 'PSDtopos' + str(cc) + '_' + raw_stem_in + '_sss_f_raw_ica_' + crop_str + fig_format)
    #     #         print('Saving figure to %s' % fname_fig)

    #     #         # Save PSD figure
    #     #         ff.savefig(fname_fig)

    #     #     # close PSD figure if in QSUB
    #     #     if close_fig:
    #     #         plt.close(ff)


    #     # fig_topo = raw.plot_psd_topo(fmin=0.5, fmax=12., dB=False, n_fft=n_fft, fig_facecolor='white')

    #     # plt.title(raw_stem_in)

    #     # fname_fig = op.join(sbj_path, 'Figures', raw_stem_in + '_sss_f_raw_psdtopo.pdf')
    #     # print('Saving figure to %s' % fname_fig)

    #     # fig_topo.savefig(fname_fig)

    #     # plt.close(fig)


    #     ### ANALYSE TARGET FREQUENCIES

    #     freq_targ = get_target_frequencies(psds_as_evo, freqs, stim_freq)

    #     # # focus on main frequencies and harmonics
    #     # n_harms = config.psd_n_harms # number of upper harmonics to consider
    #     # n_bins = config.psd_n_bins # number of neighbouring frequency bins to consider per side
    #     # n_gap = config.psd_n_gap # number of bins as "gap" between neighours (n_bins) and target frequency

    #     # # target frequencies
    #     # freq_targ = {}
    #     # freq_targ['std'] = stim_freq[0] # standard presentation frequency
    #     # freq_targ['odd'] = stim_freq[1] # oddball presentation frequency

    #     # # indices of frequencies in spectrum
    #     # freq_targ['std_idx'] = np.argmin(np.abs(freqs-freq_targ['std']))
    #     # freq_targ['odd_idx'] = np.argmin(np.abs(freqs-freq_targ['odd']))

    #     # # PSD amplitudes at target frequencies
    #     # freq_targ['std_amp'] = psds[:,freq_targ['std_idx']]
    #     # freq_targ['odd_amp'] = psds[:,freq_targ['odd_idx']]

    #     # # PSD amplitudes at neighouring bins for all channels separately
    #     # # TO DO: include harmonics

    #     # freq_targ['std_bins'] = psds[:,freq_targ['std_idx']-n_bins-n_gap:freq_targ['std_idx']-n_gap] + \
    #     #                             psds[:,freq_targ['std_idx']+n_gap+1:freq_targ['std_idx']+n_gap+n_bins+1]
    #     # freq_targ['odd_bins'] = psds[:,freq_targ['odd_idx']-n_bins-n_gap:freq_targ['odd_idx']-n_gap] + \
    #     #                             psds[:,freq_targ['odd_idx']+1+n_gap:freq_targ['odd_idx']+n_gap+n_bins+1]

    #     # # average amplitude in neighbouring bins for all channels separately
    #     # freq_targ['std_bins_avg'] = np.average(freq_targ['std_bins'], axis=1)
    #     # freq_targ['odd_bins_avg'] = np.average(freq_targ['odd_bins'], axis=1)

    #     # # standard deviation in neighbouring bins for all channels separately
    #     # freq_targ['std_bins_sd'] = np.std(freq_targ['std_bins'], axis=1)
    #     # freq_targ['odd_bins_sd'] = np.std(freq_targ['odd_bins'], axis=1)

    #     # # Z-score for target frequency vs neighbouring bins, per channel
    #     # freq_targ['std_amp_z'] = (freq_targ['std_amp'] - freq_targ['std_bins_avg'])/freq_targ['std_bins_sd']
    #     # freq_targ['odd_amp_z'] = (freq_targ['odd_amp'] - freq_targ['odd_bins_avg'])/freq_targ['odd_bins_sd']

    #     # get maximum and average values across channels for display
    #     std_z_max = freq_targ['std_amp_z'].max()
    #     idx_std = np.argmax(freq_targ['std_amp_z'])
    #     std_amp_max = freq_targ['std_amp'][idx_std]
    #     std_bins_max = freq_targ['std_bins_avg'][idx_std]
    #     std_sd_max = freq_targ['std_bins_sd'][idx_std]
    #     std_max_channel = info['ch_names'][idx_std]

    #     odd_z_max = freq_targ['odd_amp_z'].max()
    #     idx = np.argmax(freq_targ['odd_amp_z'])
    #     odd_amp_max = freq_targ['odd_amp'][idx]
    #     odd_bins_max = freq_targ['odd_bins_avg'][idx]
    #     odd_sd_max = freq_targ['odd_bins_sd'][idx]
    #     odd_max_channel = info['ch_names'][idx]

    #     print('Values for standard frequency %f (%d) in Channel %s:' % (freq_targ['std'], freq_targ['std_idx'], std_max_channel))
    #     print('Z-score %f: (%e - %e)/%e.' % (std_z_max, std_amp_max, std_bins_max, std_sd_max))

    #     print('Values for oddball frequency %f (%d) in Channel %s:' % (freq_targ['odd'], freq_targ['odd_idx'], odd_max_channel))
    #     print('Z-score %f: (%e - %e)/%e.' % (odd_z_max, odd_amp_max, odd_bins_max, odd_sd_max))

    #     # if TFR analysis requested
    #     if do_tfr:

    #         print('\nTime-frequency analysis.')

    #         # try to get large raw segment as epoch for TFR analysis
    #         epoch = mne.Epochs(raw_ori, events=events, event_id=19, tmin=config.tfr['epoch'][0],
    #                             tmax=config.tfr['epoch'][1], reject=None, preload=True)

    #         # compute TFR power without inter-trial coherence for this epoch
    #         # number of cycles is frequency-dependent
    #         n_cycles = [3 if x<7. else int(x/2.) for x in config.tfr['freqs']]            

    #         powtfr = mne.time_frequency.tfr_morlet(epoch, freqs=config.tfr['freqs'], n_cycles=n_cycles,
    #                                             return_itc=False, average=True)

    #         # baseline normalisation
    #         powtfr = powtfr.apply_baseline(baseline=(config.tfr['epoch'][0], 0.), mode='zscore')

    #         # Plot TRF results
    #         plot_TFR(powtfr, idx_std, std_max_channel, stim_freq[0], sbj_path, raw_stem_in, '', close_fig)

    #         # # get maximum value across time at standard frequency for scaling
    #         # idx_fr = np.argmin(np.abs(config.tfr['freqs']-freq_targ['std']))
    #         # vmax = powtfr.data[idx_std, idx_fr, :].max()

    #         # # plot TFR ratio for maximum channel determined above
    #         # fig = powtfr.plot(picks=std_max_channel, title=std_max_channel, vmin=-vmax, vmax=vmax)

    #         # # yticks at integer frequencies, labels every 5Hz
    #         # ticks = np.arange(int(config.tfr['freqs'][0]), int(config.tfr['freqs'][-1]), 1)
    #         # lab_freqs = np.arange(int(config.tfr['freqs'][0]), int(config.tfr['freqs'][-1]), 5)
    #         # labels = [x if x in lab_freqs else '' for x in ticks]
    #         # plt.yticks(ticks=ticks, labels=labels)

    #         # # name depends on crop_time in config.py
    #         # fname_fig = op.join(sbj_path, 'Figures', 'TFR_' + raw_stem_in + '_sss_f_raw_ica_' + crop_str + '.pdf')
    #         # print('Saving figure to %s' % fname_fig)

    #         # # Save PSD figure
    #         # fig.savefig(fname_fig)

    #         # collect results across runs
    #         powtfr_all[cond_lbl].append(powtfr)

    #     else: # if no TFR requested

    #         print('\nNot doing time-frequency analysis.')

    #         powtfr_all[cond_lbl] = []

    #     # collect results across runs
    #     freq_targ_all[cond_lbl].append(freq_targ)
    #     psds_as_evo_all[cond_lbl].append(psds_as_evo)

    # # GRAND-AVERAGE results across relevant runs
    # psds_as_evo_gm = {}
    # powtfr_gm = {}

    # for cond_lbl in freq_targ_all.keys():

    #     # average across runs
    #     if len(psds_as_evo_all[cond_lbl])>1:

    #         psds_as_evo_gm[cond_lbl] = mne.grand_average(psds_as_evo_all[cond_lbl], interpolate_bads=True)

    #     else:

    #         psds_as_evo_gm[cond_lbl] = psds_as_evo_all[cond_lbl]

    #     if do_tfr and len(powtfr_all)>1:

    #         powtfr_gm[cond_lbl] = mne.grand_average(powtfr_all[cond_lbl], interpolate_bads=True)

    #     else:

    #         powtfr_gm[cond_lbl] = powtfr_all[cond_lbl]

    #     # Plot GM results
    #     txt_label = 'GM_' + cond_lbl + '_'
    #     plot_psd_as_evo(psds_as_evo_gm[cond_lbl], sbj_path, '', txt_label, close_fig)

    #     # find the right frequencies for this condition
    #     stim_freq = freq_targ_all[cond_lbl][0]['stim_freq']

    #     freq_targ = get_target_frequencies(psds_as_evo_gm[cond_lbl], freqs, stim_freq)        

    #     # get maximum and average values across channels for display
    #     std_z_max = freq_targ['std_amp_z'].max()
    #     idx_std = np.argmax(freq_targ['std_amp_z'])
    #     std_amp_max = freq_targ['std_amp'][idx_std]
    #     std_bins_max = freq_targ['std_bins_avg'][idx_std]
    #     std_sd_max = freq_targ['std_bins_sd'][idx_std]
    #     std_max_channel = info['ch_names'][idx_std]

    #     odd_z_max = freq_targ['odd_amp_z'].max()
    #     idx = np.argmax(freq_targ['odd_amp_z'])
    #     odd_amp_max = freq_targ['odd_amp'][idx]
    #     odd_bins_max = freq_targ['odd_bins_avg'][idx]
    #     odd_sd_max = freq_targ['odd_bins_sd'][idx]
    #     odd_max_channel = info['ch_names'][idx]

    #     print('\n###\n')
    #     print('Combined runs:\n')
    #     print('Values for standard frequency %f (%d) in Channel %s:' % (freq_targ['std'], freq_targ['std_idx'], std_max_channel))
    #     print('Z-score %f: (%e - %e)/%e.' % (std_z_max, std_amp_max, std_bins_max, std_sd_max))

    #     print('Values for oddball frequency %f (%d) in Channel %s:' % (freq_targ['odd'], freq_targ['odd_idx'], odd_max_channel))
    #     print('Z-score %f: (%e - %e)/%e.' % (odd_z_max, odd_amp_max, odd_bins_max, odd_sd_max))

    #     if do_tfr:

    #         plot_TFR(powtfr_gm[cond_lbl], idx_std, std_max_channel, stim_freq[0], sbj_path, raw_stem_in, txt_label, close_fig)


    # return psds_as_evo_all, psds_as_evo_gm, powtfr_all, powtfr_gm, freq_targ_all # psds, psds_as_evo, freqs # return last raw data
    return


def psds_for_harmonics(psds, freqs, basefreq, oddfreq, n_harms, n_bins,
                       method='sum'):
    """Combine amplitudes across harmonics of oddball frequency, ignoring
    the base frequency.
    Parameters:
        psds: array
            The PSD, shape (n_channels, n_freqs)
        freqs: array-like
            The frequencies corresponding to n_freqs columns of psds.
        basefreq: float
            The base frequency (Hz).
        oddfreq: float
            The oddball frequency (Hz).
        n_harms: int
            The number of harmonics to combine.
            This does not include base frequency and its harmonics.
        n_bins: int
            Number of bins neighbouring harmonics to take into account.
        method: str
            'sum' or 'avg'
            Whether to sum or average amplitudes across harmonics.
            Defaults to 'sum'.
    Returns:
        psd_harm: array, dimension (2 * n_bins + 1)
            The combined PSDs around harmonics.
    """
    # get harmonics of oddfreq that do not overlap with harmonics of
    # basefreq
    harm_freqs = _get_valid_harmonics(freqs, basefreq, oddfreq, n_harms)

    # find indices corresponding to valid harmonic frequencies
    harm_idx = [np.abs(ff - freqs).argmin() for ff in harm_freqs]

    psd_harms = np.zeros([psds.shape[0], 2 * n_bins + 1])

    # Sum up PSDs around harmonics
    for ii in harm_idx:

        idx = np.arange(ii - n_bins, ii + n_bins + 1)

        # get PSD for bin around harmonic
        psd_now = psds[:, idx]

        psd_harms = psd_harms + psd_now

    # average if requested
    if method == 'avg':

        psd_harms = psd_harms / harm_freqs.size

    return psd_harms



def combine_harmonics(psds, freqs, basefreq, oddfreq, n_harms, method='sum'):
    """Combine amplitudes across harmonics of oddball frequency, ignoring
    the base frequency.
    Parameters:
        psds: array
            The PSD, shape (n_channels, n_freqs)
        freqs: array-like
            The frequencies corresponding to n_freqs columns of psds.
        basefreq: float
            The base frequency (Hz).
        oddfreq: float
            The oddball frequency (Hz).
        n_harms: int
            The number of harmonics to combine.
            This does not include base frequency and its harmonics.
        method: str
            'sum' or 'avg'
            Whether to sum or average amplitudes across harmonics.
            Defaults to 'sum'.
    Returns:
        sum_harm: float
            The amplitude combined across harmonics of oddball frequency.
            Also includes "zero harmonic", i.e. the oddball frequency itself.

    """
    # get harmonics of oddfreq that do not overlap with harmonics of
    # basefreq
    harm_freqs = _get_valid_harmonics(freqs, basefreq, oddfreq, n_harms)

    # find indices corresponding to valid harmonic frequencies
    harm_idx = [np.abs(ff - freqs).argmin() for ff in harm_freqs]

    print('Frequencies to be combined:')
    print(harm_freqs)

    # take only PSD values at harmonics of oddball frequency
    psds_harm = psds[:, harm_idx]

    # Average amplitudes across valid harmonic frequencies
    sum_harm = np.sum(psds_harm, axis=1)

    # if average across harmonics requested
    if method == 'avg':

        sum_harm = sum_harm / n_harms

    return sum_harm


def _get_valid_harmonics(freqs, basefreq, oddfreq, n_harms):
    """Compute harmonics of oddfreq while ignoring harmonics of basefreq.
    Parameters:
        freqs: array-like
            The frequencies corresponding to n_freqs columns of psds.
        basefreq: float
            The base frequency (Hz).
        oddfreq: float
            The oddball frequency (Hz).
        n_harms: int
            The number of harmonics to combine.
            This does not include base frequency and its harmonics.
    Returns:
        harm_freqs: array
            Valid harmonics of oddfreq that do not overlap with harmonics of
            basefreq.
    """
    # start with this list of possible harmonics of oddball
    # only up to maximum frequency in PSD
    max_freq = np.min([(n_harms + 1) * oddfreq, freqs[-1]])
    harm_freqs = np.arange(oddfreq, max_freq, oddfreq)
    # rounding necessary to make comparisons in float accurate
    harm_freqs = np.round(harm_freqs, 6)

    # base frequency and its harmonics should not be included
    check_freqs = np.arange(basefreq, (n_harms + 1) * basefreq, basefreq)
    check_freqs = np.round(check_freqs, 6)

    # TO DO: check precision of these frequencies, especially for faces!!!

    # indices of elements to be removed from harm_freqs
    del_idx = [i for i in np.arange(0, len(harm_freqs)) if harm_freqs[i] in
               check_freqs]

    print(del_idx)

    # flip because elements will be deleted, changing order for del_idx
    for idx in np.flip(del_idx):

        # delete one harmonic of oodball that coincides with base frequency
        harm_freqs = np.delete(harm_freqs, idx)

        # append new frequency
        add_freq = harm_freqs[-1] + oddfreq

        print('1 %s' % add_freq)

        # no while required because there can't be two harmonics next to each
        # other
        if add_freq in check_freqs:
            add_freq = harm_freqs[-1] + 2. * oddfreq
            print('2 %s' % add_freq)

        harm_freqs = np.append(harm_freqs, add_freq)

    return harm_freqs


def get_target_frequencies(psds_as_evo, freqs, stim_freq):
    """# Analyse PSD at target frequencies."""
# psds_as_evo: Power spectral density as Evoked object
# freqs: list of frequencies from PSD
# stim_freq: 2-item list with presentation and oddball frequencies
# returns: freq_targ, dict with results

# focus on main frequencies and harmonics
    # number of upper harmonics to consider
    n_harms = config.psd_n_harms
    # number of neighbouring frequency bins to consider per side
    n_bins = config.psd_n_bins
    # number of bins as "gap" between neighours (n_bins) and target frequency
    n_gap = config.psd_n_gap

    psds = psds_as_evo.data

    # target frequencies
    freq_targ = {}
    freq_targ['stim_freq'] = stim_freq
    freq_targ['std'] = stim_freq[0]  # standard presentation frequency
    freq_targ['odd'] = stim_freq[1]  # oddball presentation frequency

    # indices of frequencies in spectrum
    freq_targ['std_idx'] = np.argmin(np.abs(freqs - freq_targ['std']))
    freq_targ['odd_idx'] = np.argmin(np.abs(freqs - freq_targ['odd']))

    # PSD amplitudes at target frequencies
    freq_targ['std_amp'] = psds[:, freq_targ['std_idx']]
    freq_targ['odd_amp'] = psds[:, freq_targ['odd_idx']]

    # PSD amplitudes at neighouring bins for all channels separately
    # TO DO: include harmonics

    freq_targ['std_bins'] = psds[:, freq_targ['std_idx'] - n_bins
                                   - n_gap:freq_targ['std_idx'] - n_gap] + psds[:, freq_targ['std_idx'] + n_gap + 1:freq_targ['std_idx'] + n_gap + n_bins + 1]
    freq_targ['odd_bins'] = psds[:, freq_targ['odd_idx'] - n_bins
                                   - n_gap:freq_targ['odd_idx'] - n_gap] + psds[:, freq_targ['odd_idx'] + 1 + n_gap:freq_targ['odd_idx'] + n_gap + n_bins + 1]

    # average amplitude in neighbouring bins for all channels separately
    freq_targ['std_bins_avg'] = np.average(freq_targ['std_bins'], axis=1)
    freq_targ['odd_bins_avg'] = np.average(freq_targ['odd_bins'], axis=1)

    # standard deviation in neighbouring bins for all channels separately
    freq_targ['std_bins_sd'] = np.std(freq_targ['std_bins'], axis=1)
    freq_targ['odd_bins_sd'] = np.std(freq_targ['odd_bins'], axis=1)

    # Z-score for target frequency vs neighbouring bins, per channel
    freq_targ['std_amp_z'] = (freq_targ['std_amp'] - freq_targ['std_bins_avg'])/freq_targ['std_bins_sd']
    freq_targ['odd_amp_z'] = (freq_targ['odd_amp'] - freq_targ['odd_bins_avg'])/freq_targ['odd_bins_sd']

    return freq_targ


def psd_convert_to_snr(psds, n_bins, n_gap=0):
    """Compute PSD (SD) SNR with respect to neighbouring frequency bins.

    Parameters:
    psds: array
        The PSD, shape (n_channels, n_freqs).
    n_bins: int
        Number of neighbouring frequency bins to use to compute standard
        deviation. Bins will be taken from both sides at each frequency.
    n_gap: int
        Gap between target frequency and neighbouring bins.
    Returns:
    psds_snr: array
        The PSD as SNRs, shape (n_channels, n_freqs).
        SNR is amplitude divided by standard deviation.
        If z-score required, subtract mean separately using baseline
        correction.
    """
    # TO DO: faster with matrix computations?
    psds_snr = np.zeros(psds.shape)  # initialise output array
    for ff in np.arange(0, psds.shape[1]):  # for frequencies

        # take care of edges in PSD
        m = np.max([0, ff - n_bins - n_gap])
        n = np.min([psds.shape[1], ff + 1 + n_bins + n_gap])

        for cc in np.arange(0, psds.shape[0]):  # for channels

            # neighbouring elements before and after this frequency
            base_idx = np.r_[np.arange(m, ff - n_gap), np.arange(ff + 1 +
                                                                 n_gap, n)]
            # get baseline amplitudes
            baseline = psds[cc, base_idx]

            # Compute SNR for one frequency
            psds_snr[cc, ff] = psds[cc, ff] / np.std(baseline)

    return psds_snr


def psd_correct_baseline(psds, n_bins, n_gap=0):
    """BASELINE-CORRECT PSDs with neighbouring frequency bins.
    Parameters:
    psds: array
        PSD, shape (n_channels x n_freqs), from mne.time_frequency.psd_welch
    n_gap: int
        Gap between target frequency and neighbouring bins.
    Returns:
    psds_base: array
        Baseline-corrected PSD, shape (n_channels x n_freqs)
    """
    psds_base = np.zeros(psds.shape)  # initialise output array
    for ff in np.arange(0, psds.shape[1]):  # for frequencies

        # take care of edges in PSD
        m = np.max([0, ff - n_bins - n_gap])
        n = np.min([psds.shape[1], ff + n_bins + n_gap])

        for cc in np.arange(0, psds.shape[0]):  # for channels

            # neighbouring elements before and after this frequency
            baseline = np.r_[psds[cc, m:ff - n_gap], psds[cc, ff + n_gap + 1:n + 1]]

            # baseline-correct at this frequency
            psds_base[cc, ff] = psds[cc, ff] - np.average(baseline)

    return psds_base


# Plot PSDs in Evoked format
def plot_psd_as_evo(psds_as_evo, sbj_path, txt_label='', close_fig=1,
                    scalings=dict(eeg=1e6, grad=1e13, mag=1e15)):
    """Plot PSD disguised as instance of Evoked.
    Parameters:
        psds_as_evo: The PSD as Evoked object
        sbj_path: path where "Figures" sub-directory is
        txt_label: string to make filename specific (at beginning)
        close_fig: close figure (1) or not (0)
        scalings: dict, scalings for eeg/grad/mag
    Returns:
        figs: list
        The list of pyplot figures created.
    """
    # keep a copy for scaling below
    psds_tmp = deepcopy(psds_as_evo)

    # CROP PSD for display
    psds_as_evo.crop(tmin=config.crop_times[0], tmax=config.crop_times[1])

    # plot topographic maps for different frequencies
    data = psds_as_evo.data

    # EEG present?
    is_eeg = psds_as_evo.__contains__('eeg')

    # quick hack to get scaling approximately right
    # avoid first 1Hz in scaling
    psds_tmp.crop(config.crop_times[0], tmax=config.crop_times[1])

    # Default mne-python scalings for evo.plot(), just for clarity
    # scalings = dict(eeg=1e6, grad=1e13, mag=1e15)
    # show y-axis with original values
    # scalings = dict(eeg=1., grad=1., mag=1.)

    # units for y-axis of PSD plots
    units = {'mag': r'amp/$\sqrt{Hz}$', 'grad': r'amp/$\sqrt{Hz}$'}

    ch_types = ['mag', 'grad']

    eeg, meg = False, True
    if is_eeg:  # only if EEG in fiff-file
        eeg = True
        units['eeg'] = r'amp/$\sqrt{Hz}$'
        ch_types.append('eeg')

    # plot y-axis range (can be negative after baseline-correction)
    ylim = {'mag': [np.min(psds_tmp.data[2:306:3, :]) * scalings['mag'],
                    np.max(psds_tmp.data[2:306:3, :]) * scalings['mag']],
            'grad': [np.min([psds_tmp.data[0:306:3, :],
                            psds_tmp.data[1:306:3, :]]) * scalings['grad'],
                     np.max([psds_tmp.data[0:306:3, :],
                            psds_tmp.data[1:306:3, :]]) * scalings['grad']]}
    if is_eeg:
        ylim['eeg'] = [np.min(psds_tmp.data[306:376,:] * scalings['eeg']),
                       np.max(psds_tmp.data[306:376,:] * scalings['eeg'])]

    print('Upper limits for Mag: %e, Grad: %e ' %
          (ylim['mag'][1], ylim['grad'][1]))

    if is_eeg:
        print('EEG: %e.\n' % ylim['eeg'][1])

    # Different target frequency for faces
    if 'face' in txt_label:

        ftimes = config.topo_times['faces']

    else:

        ftimes = config.topo_times['words']

    ts_args = dict(spatial_colors=True, scalings=scalings, units=units,
                   ylim=ylim, time_unit='ms')

    figs = []
    for ch_type in ch_types:

        topomap_args = dict(scalings=scalings, time_format='%.2f Hz',
                            time_unit='ms')

        fig = psds_as_evo.plot_joint(times=ftimes, title=txt_label,
                                     ts_args=ts_args, picks=ch_type,
                                     topomap_args=topomap_args)

        figs.append(fig)

    for (cc, ff) in zip(ch_types, figs):

        for fig_format in ['.pdf']:  # jpg doesn't work, png does

            fname_fig = op.join(sbj_path, 'Figures', txt_label + str(cc) +
                                fig_format)
            print('Saving figure to %s' % fname_fig)

            # Save PSD figure
            ff.savefig(fname_fig)

        # close PSD figure if in QSUB
        if close_fig:
            plt.close(ff)

    return figs


# plot TRF results
def plot_TFR(powtfr, idx_std, std_max_channel, freq_std, sbj_path, raw_stem_in, txt_label, close_fig):
# powtfr: TFRAverage object
# idx_std: index to standard frequency
# std_max_channel: name of channel with maximum amplitude at standard frequency
# freq_std: standard presentation frequency
# sbj_path: path where "Figures" sub-directory is
# raw_stem_in: part of the filename, for PDF filename
# txt_label: string to make filename more specific
# close_fig: close figure or not

    # get maximum value across time at standard frequency for scaling
    idx_fr = np.argmin(np.abs(config.tfr['freqs']-freq_std))
    vmax = powtfr.data[idx_std, idx_fr, :].max()

    # plot TFR ratio for maximum channel determined above
    fig = powtfr.plot(picks=std_max_channel, title=std_max_channel, vmin=-vmax, vmax=vmax)

    # yticks at integer frequencies, labels every 5Hz
    ticks = np.arange(int(config.tfr['freqs'][0]), int(config.tfr['freqs'][-1]), 1)
    lab_freqs = np.arange(int(config.tfr['freqs'][0]), int(config.tfr['freqs'][-1]), 5)
    labels = [x if x in lab_freqs else '' for x in ticks]
    plt.yticks(ticks=ticks, labels=labels)

    crop_str = str(int(10000.*config.crop_times[0])) + '_' + str(int(10000.*config.crop_times[1]))

    # name depends on crop_time in config.py
    fname_fig = op.join(sbj_path, 'Figures', 'TFR_' + txt_label + raw_stem_in + '_sss_f_raw_ica_' + crop_str + '.pdf')
    print('Saving figure to %s' % fname_fig)

    # Save PSD figure
    fig.savefig(fname_fig)

    # close PSD figure if in QSUB
    if close_fig:
        plt.close(fig)


# get all input arguments except first
if len(sys.argv)==1:

    sbj_ids = np.arange(0,len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    # raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
    run_PSD_raw(ss)

print('Done.')
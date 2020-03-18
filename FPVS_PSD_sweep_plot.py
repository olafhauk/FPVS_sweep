#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
Plot results of FPVS Frequency Sweep, produced by FPVS_PSD_sweep_compute.py.

Average raw data from FPVS_get_sweeps.py.
Compute z-scores.
Compute TFR if specified.
==========================================

OH, October 2019
"""

### NOT DONE YETS

import sys

from os import path as op
import numpy as np

# import matplotlib
# matplotlib.use('Agg') #  for running graphics on cluster ### EDIT

from matplotlib import pyplot as plt

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

print('Sunshine.')


def run_PSD_raw(sbj_id):
    """Compute spectra for one subject."""
    # initialise html report for one subject
    report = Report(subject=str(sbj_id), title='FPVS PSDs')

    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # path to sub-directory for figures
    figs_path = op.join(sbj_path, figs_dir)

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    # get condition names and frequency names
    conds = []  # names of conditions
    for raw_stem_in in sss_map_fname[1][2:]:

        conds.append(raw_stem_in[:4])

    conds = np.unique(conds)

    # Base frequencies for frequency sweep for words (not faces)
    freqs_all = [str(ff) for ff in config.fpvs_freqs]

    print('Frequencies used: ')
    print(freqs_all)

    # initialise sum across harmonics for conditions
    sum_harms = {}
    for cond in conds:

        sum_harms[cond] = {}

    # Go through conditions and frequencies
    for cond in conds:  # conditions

        # read Evoked objects for all frequencies per condition
        print('Reading results from Evoked files.')

        # separate filename prefixes for ICAed and non-ICAed data
        prefix = ''
        if 'ica' in config.raw_ICA_suff:
            prefix = 'ICA'

        fname_evo = op.join(sbj_path, '%sPSDTopo_%s%s' % (prefix, cond,
                                                          '-ave.fif'))
        print(fname_evo)
        psds_as_evo = mne.read_evokeds(fname_evo)

        fname_evo = op.join(sbj_path, '%sPSDTopoZ_%s%s' % (prefix, cond,
                                                           '-ave.fif'))
        print(fname_evo)
        psds_z_as_evo = mne.read_evokeds(fname_evo)

        fname_evo = op.join(sbj_path, '%sPSDHarm_%s%s' % (prefix, cond,
                                                          '-ave.fif'))
        print(fname_evo)
        psd_harm_as_evo = mne.read_evokeds(fname_evo)

        fname_evo = op.join(sbj_path, '%sPSDHarmBase_%s%s' % (prefix, cond,
                                                              '-ave.fif'))
        print(fname_evo)
        psd_harm_base_as_evo = mne.read_evokeds(fname_evo)

        fname_evo = op.join(sbj_path, '%sPSDSumTopoOdd_%s%s' % (prefix, cond,
                                                                '-ave.fif'))
        print(fname_evo)
        sum_odd_as_evo = mne.read_evokeds(fname_evo)

        fname_evo = op.join(sbj_path, '%sPSDSumTopoBase_%s%s' % (prefix, cond,
                                                                 '-ave.fif'))
        print(fname_evo)
        sum_base_as_evo = mne.read_evokeds(fname_evo)

        # Establish channel types present in these data
        chtypes = ['mag', 'grad', 'eeg']  # all possible channel types

        for [ci, ch_type] in enumerate(chtypes):

            if not psds_as_evo[0].__contains__(ch_type):

                del(chtypes[chtypes.index[ch_type]])

        if cond == 'face':  # no frequency sweep for faces

            freqs = ['6.0']  # base frequency for this condition (Hz as string)

            freq_odd = 1.2  # oddball frequency for this condition (Hz)

        else:  # for all word condition, use all sweep frequencies

            # base frequencies for this condition (Hz as string)
            freqs = freqs_all

            freq_odd = 1.0  # oddball frequency the same for all sweeps

        for [fi, freq] in enumerate(freqs):  # frequencies

            basefreq = float(freq)  # hack, float-to-string-to-float-again

            # label for condition and base frequency
            label_str = '%s_%s' % (cond, ''.join(freq.split('.')))

            # Plot PSD as spectrum plus topographies (plot_joint())
            print('Plotting PSDs for %s.' % label_str)

            # Plot z-scored PSD
            evoked = psds_z_as_evo[fi]

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

            file_label = 'PSDTopoZ_%s' % label_str
            figs = Ff.plot_psd_as_evo(evoked, sbj_path, picks=None,
                                      txt_label=file_label,
                                      close_fig=close_fig,
                                      scalings=unit_scalings)

            for [fig, chtype] in zip(figs, chtypes):

                sec_label = '%s_%s_Z_' % (evoked.comment, chtype)

                report.add_figs_to_section(fig, sec_label, section=sec_label,
                                           scale=1)

            # Plot for peak channels (base frequency) without topographies
            fig = evoked.plot(spatial_colors=True, picks=peak_ch_names_base,
                              scalings=unit_scalings, gfp=True, time_unit='s')

            fname_fig = op.join(figs_path, 'PSDTopoZPeakbase_%s.pdf' % label_str)

            fig.savefig(fname_fig)

            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            # Plot for peak channels (oddball frequency) without topographies
            fig = evoked.plot(spatial_colors=True, picks=peak_ch_names_odd,
                              scalings=unit_scalings, gfp=True, time_unit='s')

            fname_fig = op.join(figs_path, 'PSDTopoZPeakodd_%s.pdf' % label_str)

            fig.savefig(fname_fig)

            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            plt.close('all')

            # Plot PSD across harmonics for oddball frequency
            evoked = psd_harm_as_evo[fi]

            # Plotting PSDs across harmonics (base frequency)
            fig = evoked.plot(spatial_colors=True, picks=peak_ch_names_base,
                              scalings=unit_scalings, gfp=True, time_unit='s')

            fname_fig = op.join(figs_path, 'PSDHarmOddPeakbase_%s_%s.pdf'
                                % (cond, freq))

            fig.savefig(fname_fig)

            sec_label = evoked.comment
            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            # Plotting PSDs across harmonics (base frequency)
            fig = evoked.plot(spatial_colors=True, picks=peak_ch_names_odd,
                              scalings=unit_scalings, gfp=True, time_unit='s')

            fname_fig = op.join(figs_path, 'PSDHarmOddPeakodd_%s_%s.pdf'
                                % (cond, freq))

            fig.savefig(fname_fig)

            sec_label = evoked.comment
            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            plt.close('all')

            # Plot PSD
            evoked = psds_as_evo[fi]

            file_label = 'PSDTopo_%s' % label_str

            figs = Ff.plot_psd_as_evo(evoked, sbj_path, picks=None,
                                      txt_label=file_label,
                                      close_fig=close_fig)

            for [fig, chtype] in zip(figs, chtypes):

                sec_label = '%s_%s_' % (evoked.comment, chtype)

                report.add_figs_to_section(fig, sec_label, section=sec_label,
                                           scale=1)

            # Plot for peak channels (base frequency) without topographies
            fig = evoked.plot(spatial_colors=True, picks=peak_ch_names_base,
                              scalings=unit_scalings, gfp=True, time_unit='s')

            fname_fig = op.join(figs_path, 'PSDTopoPeakbase_%s.pdf' % label_str)

            fig.savefig(fname_fig)

            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            # Plot for peak channels (oddball frequency) without topographies
            fig = evoked.plot(spatial_colors=True, picks=peak_ch_names_odd,
                              scalings=unit_scalings, gfp=True, time_unit='s')

            fname_fig = op.join(figs_path, 'PSDTopoPeakodd_%s.pdf' % label_str)

            fig.savefig(fname_fig)

            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            plt.close(fig)

            # Plot PSD across harmonics for peak channels for base frequency
            evoked = psd_harm_base_as_evo[fi]

            # Plotting PSDs across harmonics
            fig = evoked.plot(spatial_colors=True, picks=peak_ch_names_base,
                              scalings=unit_scalings, gfp=True, time_unit='s')

            fname_fig = op.join(figs_path, 'PSDHarmBasePeakbase_%s_%s.pdf'
                                % (cond, freq))

            fig.savefig(fname_fig)

            sec_label = evoked.comment
            report.add_figs_to_section(fig, sec_label, section=sec_label,
                                       scale=1)

            plt.close('all')

            # Plot PSD topography across harmonics for oddball frequency
            evoked = sum_odd_as_evo[fi]

            # Note: also for oddball frequency the "latency" is the base
            # frequency, because that's our experimental manipulation
            times = [basefreq]

            # Filename stem for figure; channel type to be added later
            fname_fig = op.join(figs_path, 'PSDSumTopoOdd_%s_%s' % (cond, freq))

            # For html report section label
            sec_label = evoked.comment

            figs = Ff.plot_evo_topomap(evoked, times, chtypes, fname_fig)

            for [fig, chtype] in zip(figs, chtypes):

                    sec_label = evoked.comment

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

            # Plot PSD topography across harmonics for base frequency
            evoked = sum_base_as_evo[fi]

            times = [basefreq]

            # Filename stem for figure; channel type to be added later
            fname_fig = op.join(figs_path, 'PSDSumTopoBase_%s_%s' % (cond, freq))

            # For html report section label
            sec_label = evoked.comment

            figs = Ff.plot_evo_topomap(evoked, times, chtypes, fname_fig)

            for [fig, chtype] in zip(figs, chtypes):

                    sec_label = evoked.comment

                    report.add_figs_to_section(fig, sec_label,
                                               section=sec_label, scale=1)

        # Save HTML report
        fname_report = op.join(figs_path, str(sbj_id) +
                               '_report.html')

        report.save(fname_report, overwrite=True, open_browser=False)

        # In case someone was forgotten
        plt.close('all')

    return


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    # raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
    run_PSD_raw(ss)

# Functions for FPVS EEG/MEG analysis
# OH

import os
from os import path as op
import numpy as np

from copy import deepcopy

from matplotlib import pyplot as plt

from importlib import reload

# from mne import EvokedArray
from mne.evoked import EvokedArray
from mne.source_estimate import SourceEstimate

# FPVS-specific parameters
import config_sweep as config
reload(config)

# for some plots of SNRs
unit_scalings = dict(eeg=1., mag=1., grad=1.)


def peak_channels_evoked(evoked, peak_times, ch_types=None, n_chan=1):
    """Reduce evoked data to peak channels per channel type.

    Parameters:
    evoked: instance of Evoked
        The evoked data for which to find peak channels.
    peak_times: list
        The latencies (s) at which to find peak channels.
    ch_types: list of string
        Channel types to be considered. 'mag' | 'grad' | 'eeg'.
        If None use all channel types in evoked.
    n_chan: int
        The number of peak channels to return per channel type.
        Default: 1.

    Returns:
    peak_ch_names: list of dict of list of strings
        The list of names of peak channels per channel type.
        [peak_times][ch_types][names]
    """
    # all possible channel types
    ch_types = ['mag', 'grad', 'eeg']

    for [ci, ch_type] in enumerate(ch_types):

        if not evoked.__contains__(ch_type):

            del(ch_types[ci])

    # indices to specified peak latencies
    peak_indices = evoked.time_as_index(peak_times)

    peak_ch_names = []

    for peak_idx in peak_indices:

        peak_ch_names.append({})

        for ch_type in ch_types:

            peak_ch_names[-1][ch_type] = []

            # channel types will be dropped
            evo_copy = deepcopy(evoked)

            if ch_type in ['mag', 'grad']:

                evo_chtype = evo_copy.pick_types(meg=ch_type, eeg=False)

            elif ch_type == 'eeg':

                evo_chtype = evo_copy.pick_types(meg=False, eeg=True)

            else:

                print('Channel type ''%s'' not recognised.' % ch_type)
                return

            topodata = np.abs(evo_chtype.data[:, peak_idx])

            # get channel indices sorted by amplitudes (ascending order)
            sort_idx = np.argsort(topodata)

            # indices to maximum channels
            max_indices = sort_idx[-1:-n_chan - 1:-1]

            # get peak channel names
            chs_peak = [evo_chtype.ch_names[i] for i in max_indices]

            peak_ch_names[-1][ch_type] = chs_peak

    return peak_ch_names


def plot_evo_topomap(evoked, times, chtypes, fname_fig):
    """Plot topographies for sum across harmonics.

    Parameters:
    evoked: Evoked instance
        The topography to plot.
    times: list
        List with one number for time axis.
    chtypes: list of str
        The channel types to plot.
    fname_fig: str
        Filename stem for figure. Channel type will be added.

    Returns:
    figs: list
        List of plot_topomap figures.
    """
    figs = []  # collect figures

    for chtype in chtypes:

        print(chtype)

        # plot topography for one channel type at a time
        fig = evoked.plot_topomap(times=times, ch_type=chtype, vmin=0.,
                                  time_format='',
                                  scalings=unit_scalings[chtype],
                                  units='Z', show=False)

        # filename for figure
        fname_fig_ch = '%s_%s.jpg' % (fname_fig, chtype)

        print('Saving figure for combined topography across harmonics to: %s.'
              % fname_fig_ch)

        fig.savefig(fname_fig_ch)

        figs.append(fig)

    return figs


def psds_across_harmonics(psd, freqs, basefreq, oddfreq, n_harms, n_bins,
                          n_gap=0, method='sum'):
    """Combine across harmonics of oddball frequency, without base frequency.

    Parameters:
        psd: instance of Evoked or Source Estimate
            The PSD, data of shape (n_channels, n_freqs)
        freqs: array-like
            The frequencies corresponding to n_freqs columns of psd.
        basefreq: float | None
            The base frequency (Hz). Its harmonics are to be excluded.
            If none, use all harmonics of oddfreq.
        oddfreq: float
            The oddball frequency (Hz).
        n_harms: int
            The number of harmonics to combine.
            This does not include base frequency and its harmonics.
        n_bins: int
            Number of bins neighbouring harmonics to take into account.
        n_gap: int
            Gap between target frequency and neighbouring bins.
        method: str
            'sum' or 'avg'
            Whether to sum or average amplitudes across harmonics.
            Defaults to 'sum'.
    Returns:
        psd_harms: instance of Evoked or Source Estimate,
            data of dimension (2 * n_bins + 1).
            The combined PSDs around harmonics.
        topo: instance of Evoked or Source Estimate
            Topography of the combined responses across harmonics at centre
            frequency (0 Hz in psd_harms).
        topos: instance of Evoked or Source Estimate
            Topographies across harmonics.
        freqs_harm: list of float
            Frequencies of harmonics taken into account.

    """
    # get data into numpy array, works for Evoked and SourceEstimate
    data = psd.data

    if basefreq is not None:
        # get harmonics of oddfreq that do not overlap with harmonics of
        # basefreq
        freqs_harm = _get_valid_harmonics(freqs, basefreq, oddfreq, n_harms)

    else:

        freqs_harm = np.arange(oddfreq, (n_harms + 1) * oddfreq, oddfreq)

    # find indices corresponding to valid harmonic frequencies
    harm_idx = [np.abs(ff - freqs).argmin() for ff in freqs_harm]

    # initialise sum of PSD-segments across harmonics
    data_harms = np.zeros([data.shape[0], 2 * n_bins + 2 * n_gap + 1])

    # initialise topographies across harmonics
    topos_mat = np.zeros([data.shape[0], len(freqs_harm)])

    # Sum up PSDs around harmonics
    for (ii, iii) in enumerate(harm_idx):

        idx = np.arange(iii - n_bins - n_gap, iii + n_bins + n_gap + 1)

        # get PSD for bin around harmonic
        data_now = data[:, idx]

        data_harms = data_harms + data_now

        # collect topography at harmonic
        topos_mat[:, ii] = data_now[:, n_bins + n_gap]

    # average if requested
    if method == 'avg':

        data_harms = data_harms / freqs_harm.size

    freq_resol = freqs[1] - freqs[0]

    tmin = -(n_bins + n_gap) * freq_resol  # include baseline

    # put processed data into Evoked or SourceEstimate,
    # depending on input
    if type(psd) is SourceEstimate:

        vertices = [psd.lh_vertno, psd.rh_vertno]

        tstep = freq_resol

        psd_harms = SourceEstimate(data=data_harms, vertices=vertices,
                                   tmin=tmin, tstep=tstep)

        # topography for combined responses at centre frequency
        t = data_harms[:, n_bins + n_gap]

        topo = SourceEstimate(
            data=t[np.newaxis, :].T, vertices=vertices,
            tmin=0., tstep=0.001)

        # topographies for all harmonics
        topos = SourceEstimate(
            data=topos_mat, vertices=vertices, tmin=0., tstep=0.001)

    elif type(psd) is EvokedArray:

        info = psd.info

        info['sfreq'] = 1. / freq_resol  # to display samples as time points

        nave = psd.nave

        psd_harms = EvokedArray(data_harms, info, tmin=tmin, nave=nave)

        # topography for combined responses at centre frequency
        t = data_harms[:, n_bins + n_gap]

        print(t.shape)

        topo = EvokedArray(
            t[np.newaxis, :].T, info, tmin=0., nave=nave)

        # topographies for all harmonics
        topos = EvokedArray(
            topos_mat, info, tmin=0.)

    else:

        print('Type of ''psd'' not known (%s)' % type(psd))

    return psd_harms, topo, topos, freqs_harm


def combine_harmonics_topos(psd, freqs, basefreq, oddfreq, n_harms,
                            method='sum'):
    """Combine topographies across harmonics of oddball frequency.

    Parameters:
        psd: instance of Evoked or Source Estimate
            The PSD, data of shape (n_channels, n_freqs)
        freqs: array-like
            The frequencies corresponding to n_freqs columns of psd.
        basefreq: float | None
            The base frequency (Hz). Its harmonics are to be excluded.
            If none, use all harmonics of oddfreq.
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
        psd_sum_harm: instance of Evoked or Source Estimate
            The amplitude combined across harmonics of oddball frequency.
            Also includes "zero harmonic", i.e. the oddball frequency itself.
            Multiples of base frequency are not included.
        topos_harm: instance of Evoked or Source Estimate
            The topographies of harmonics that went into the sum or average.
        freqs_harm: numpy array
            The frequencies corresponding to the harmonics that went into the
            sum or average.
    """
    # get data as array, works for Evoked and SourceEstimate
    data = psd.data

    if basefreq is not None:
        # get harmonics of oddfreq that do not overlap with harmonics of
        # basefreq
        freqs_harm = _get_valid_harmonics(freqs, basefreq, oddfreq, n_harms)

    else:

        freqs_harm = np.arange(oddfreq, (n_harms + 1) * oddfreq, oddfreq)

    # find indices corresponding to valid harmonic frequencies
    harm_idx = [np.abs(ff - freqs).argmin() for ff in freqs_harm]

    print('Frequencies to be combined:')
    print(freqs_harm)

    # take only PSD values at harmonics of oddball frequency
    data_harm = data[:, harm_idx]

    # Average amplitudes across valid harmonic frequencies
    sum_harm = np.sum(data_harm, axis=1)

    # if average across harmonics requested
    if method == 'avg':

        sum_harm = sum_harm / n_harms

    # insert data into Evoked object

    # needs another dimension for Evoked or SourceEstimate object
    sum_harm = np.expand_dims(sum_harm, 1)

    # put processed data into Evoked or SourceEstimate,
    # depending on input
    if type(psd) is SourceEstimate:

        vertices = [psd.lh_vertno, psd.rh_vertno]

        tmin = 0.

        tstep = 0.001

        psd_sum_harm = SourceEstimate(data=sum_harm, vertices=vertices,
                                      tmin=tmin, tstep=tstep)

        # as Evoked for return
        topos_harm = SourceEstimate(data=data_harm, vertices=vertices,
                                    tmin=tmin, tstep=tstep)

    elif type(psd) is EvokedArray:

        nave = psd.nave

        psd_sum_harm = EvokedArray(sum_harm, psd.info, tmin=0., nave=nave)

        topos_harm = EvokedArray(data_harm, psd.info, tmin=0., nave=nave)

    else:

        print('Type of ''psd'' not known (%s)' % type(psd))

    return psd_sum_harm, topos_harm, freqs_harm


def _get_valid_harmonics(freqs, basefreq, oddfreq, n_harms):
    """Compute harmonics of oddfreq while ignoring harmonics of basefreq.

    Parameters:
        freqs: array-like
            The frequencies corresponding to n_freqs columns of psd.
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

    # indices of elements to be removed from harm_freqs
    del_idx = [i for i in np.arange(0, len(harm_freqs)) if harm_freqs[i] in
               check_freqs]

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


# NOT USED
def get_target_frequencies(psds_as_evo, freqs, stim_freq):
    """Analyse PSD at target frequencies."""
    # psds_as_evo: Power spectral density as Evoked object
    # freqs: list of frequencies from PSD
    # stim_freq: 2-item list with presentation and oddball frequencies
    # returns: freq_targ, dict with results

# focus on main frequencies and harmonics
    # number of neighbouring frequency bins to consider per side
    n_bins = config.psd_n_bins
    # number of bins as "gap" between neighours (n_bins) and target frequency
    n_gap = config.psd_n_gap

    data = psds_as_evo.data

    # target frequencies
    freq_targ = {}
    freq_targ['stim_freq'] = stim_freq
    freq_targ['std'] = stim_freq[0]  # standard presentation frequency
    freq_targ['odd'] = stim_freq[1]  # oddball presentation frequency

    # indices of frequencies in spectrum
    freq_targ['std_idx'] = np.argmin(np.abs(freqs - freq_targ['std']))
    freq_targ['odd_idx'] = np.argmin(np.abs(freqs - freq_targ['odd']))

    # PSD amplitudes at target frequencies
    freq_targ['std_amp'] = data[:, freq_targ['std_idx']]
    freq_targ['odd_amp'] = data[:, freq_targ['odd_idx']]

    # PSD amplitudes at neighouring bins for all channels separately
    # TO DO: include harmonics

    freq_targ['std_bins'] = (
        data[:, freq_targ['std_idx'] - n_bins - n_gap:freq_targ['std_idx'] -
             n_gap] +
        data[:, freq_targ['std_idx'] + n_gap + 1:freq_targ['std_idx'] +
             n_gap + n_bins + 1]
    )
    freq_targ['odd_bins'] = (
        data[:, freq_targ['odd_idx'] - n_bins - n_gap:freq_targ['odd_idx'] -
             n_gap] +
        data[:, freq_targ['odd_idx'] + 1 + n_gap:freq_targ['odd_idx'] +
             n_gap + n_bins + 1]
    )

    # average amplitude in neighbouring bins for all channels separately
    freq_targ['std_bins_avg'] = np.average(freq_targ['std_bins'], axis=1)
    freq_targ['odd_bins_avg'] = np.average(freq_targ['odd_bins'], axis=1)

    # standard deviation in neighbouring bins for all channels separately
    freq_targ['std_bins_sd'] = np.std(freq_targ['std_bins'], axis=1)
    freq_targ['odd_bins_sd'] = np.std(freq_targ['odd_bins'], axis=1)

    # Z-score for target frequency vs neighbouring bins, per channel
    freq_targ['std_amp_z'] = (
        (freq_targ['std_amp'] - freq_targ['std_bins_avg']) /
        freq_targ['std_bins_sd']
    )
    freq_targ['odd_amp_z'] = (
        (freq_targ['odd_amp'] - freq_targ['odd_bins_avg']) /
        freq_targ['odd_bins_sd']
    )

    return freq_targ


def psd_z_score(psd, n_bins, mode='z', n_gap=0, minmax=True):
    """Z-score PSD with respect to neighbouring frequency bins.

    Parameters:
    psd: instance of Evoked or SourceEstimate
        The PSD, data of dimension (n_channels, n_freqs).
    n_bins: int
        Number of neighbouring frequency bins to use to compute standard
        deviation. Bins will be taken from each side at each frequency,
        thus altogether 2*n_bins bins will be used.
    mode: str ('z' | 'baseline' | 'snr')
        Whether to do subtract baseline ('baseline'), divide by standard
        deviation ('snr'), or both ('z', first baseline then snr).
        Defaults to 'z'.
    n_gap: int
        Gap between target frequency and neighbouring bins.
    minmax: Bool
        Whether or not to remove minimum and maximum values from z-scoring.
    Returns:
    psd_z: instance of Evoked or SourceEstimate
        The transformed PSD, shape (n_channels, n_freqs).
        'z' subtracts baseline and divides by standard deviation;
        snr' is amplitude divided by standard deviation;
        'baseline' only subtracts baseline.
    """
    # get PSD as numpy array (works for Evoked and SourceEstimate)
    data = psd.data

    # initialise output
    data_z = deepcopy(data)

    # will contain baseline values to be subtracted
    base_mat = np.zeros(data.shape)

    # will contain standard deviation to be divided by
    sd_mat = np.zeros(data.shape)

    # Compute baseline values and standard deviations first as matrices,
    # then subtract/divide on whole matrices
    # hopefully faster

    for ff in np.arange(0, data.shape[1]):  # for frequencies

        # take care of edges in PSD
        m = np.max([0, ff - n_bins - n_gap])
        n = np.min([data.shape[1], ff + 1 + n_bins + n_gap])

        # neighbouring elements before and after this frequency
        base_idx = np.r_[np.arange(m, ff - n_gap), np.arange(ff + 1 +
                                                             n_gap, n)]

        for cc in np.arange(0, data.shape[0]):  # for channels or vertices

            # get baseline amplitudes
            baseline = data[cc, base_idx]

            # indices of minimum and maximum baseline value
            minidx = baseline.argmin()
            maxidx = baseline.argmax()

            # remove min/max values from baseline
            np.delete(baseline, (minidx, maxidx))

            # avoiding too many comparisons
            if mode is 'baseline':

                # average taken later, since always same number of elements
                base_mat[cc, ff] = np.mean(baseline)

            elif mode == 'snr':

                sd_mat[cc, ff] = np.std(baseline)

            elif mode == 'z':

                base_mat[cc, ff] = np.mean(baseline)

                sd_mat[cc, ff] = np.std(baseline)

    # subtract average of baseline
    if mode in ['z', 'baseline']:

        data_z = data_z - base_mat

    # divided by standard deviation of baseline
    if mode in ['z', 'snr']:

        data_z = data_z / sd_mat

    # put processed data into Evoked or SourceEstimate,
    # depending on input
    if type(psd) is SourceEstimate:

        vertices = [psd.lh_vertno, psd.rh_vertno]

        tmin = psd.times[0]

        if len(psd.times) > 1:

            tstep = psd.times[1] - psd.times[0]

        else:  # in case only one sample present

            tstep = 0.001

        psd_z = SourceEstimate(data=data_z, vertices=vertices, tmin=tmin,
                               tstep=tstep)

    elif type(psd) is EvokedArray:

        tmin = psd.times[0]

        nave = psd.nave

        psd_z = EvokedArray(data_z, psd.info, tmin=tmin, nave=nave)

    else:

        print('Type of ''psd'' not know (%s)' % type(psd))

    return psd_z


# BEFORE "optimisation" of baseline and z-score
# def psd_z_score(psd, n_bins, mode='z', n_gap=0, minmax=True):
#     """Z-score PSD with respect to neighbouring frequency bins.

#     Parameters:
#     psd: instance of Evoked or SourceEstimate
#         The PSD, data of dimension (n_channels, n_freqs).
#     n_bins: int
#         Number of neighbouring frequency bins to use to compute standard
#         deviation. Bins will be taken from each side at each frequency,
#         thus altogether 2*n_bins bins will be used.
#     mode: str ('z' | 'baseline' | 'snr')
#         Whether to do subtract baseline ('baseline'), divide by standard
#         deviation ('snr'), or both ('z', first baseline then snr).
#         Defaults to 'z'.
#     n_gap: int
#         Gap between target frequency and neighbouring bins.
#     minmax: Bool
#         Whether or not to remove minimum and maximum values from z-scoring.
#     Returns:
#     psd_z: instance of Evoked or SourceEstimate
#         The transformed PSD, shape (n_channels, n_freqs).
#         'z' subtracts baseline and divides by standard deviation;
#         snr' is amplitude divided by standard deviation;
#         'baseline' only subtracts baseline.
#     """
#     # TO DO: faster with matrix computations?
#     # get PSD as numpy array (works for Evoked and SourceEstimate)
#     data = psd.data

#     # initialise output
#     data_z = deepcopy(data)

#     for ff in np.arange(0, data.shape[1]):  # for frequencies

#         # take care of edges in PSD
#         m = np.max([0, ff - n_bins - n_gap])
#         n = np.min([data.shape[1], ff + 1 + n_bins + n_gap])

#         for cc in np.arange(0, data.shape[0]):  # for channels or vertices

#             # neighbouring elements before and after this frequency
#             base_idx = np.r_[np.arange(m, ff - n_gap), np.arange(ff + 1 +
#                                                                  n_gap, n)]
#             # get baseline amplitudes
#             baseline = data[cc, base_idx]

#             # indices of minimum and maximum baseline value
#             minidx = baseline.argmin()
#             maxidx = baseline.argmax()

#             # remove min/max values from baseline
#             np.delete(baseline, (minidx, maxidx))

#             if mode in ['z', 'baseline']:

#                 # baseline-correct at this frequency
#                 data_z[cc, ff] = data_z[cc, ff] - np.average(baseline)

#             if mode in ['z', 'snr']:

#                 # Compute SNR for one frequency
#                 data_z[cc, ff] = data_z[cc, ff] / np.std(baseline)

#     # put processed data into Evoked or SourceEstimate,
#     # depending on input
#     if type(psd) is SourceEstimate:

#         vertices = [psd.lh_vertno, psd.rh_vertno]

#         tmin = psd.times[0]

#         if len(psd.times) > 1:

#             tstep = psd.times[1] - psd.times[0]

#         else:  # in case only one sample present

#             tstep = 0.001

#         psd_z = SourceEstimate(data=data_z, vertices=vertices, tmin=tmin,
#                                tstep=tstep)

#     elif type(psd) is EvokedArray:

#         psd_z = EvokedArray(data_z, psd.info)

#     else:

#         print('Type of ''psd'' not know (%s)' % type(psd))

#     return psd_z

### copy of old version that works with numpy arrays
# # def psd_convert_to_snr(psds, n_bins, n_gap=0):
# def psd_z_score(psds, n_bins, mode='z', n_gap=0, minmax=True):
#     """Compute PSD (SD) SNR with respect to neighbouring frequency bins.

#     Parameters:
#     psds: array
#         The PSD, shape (n_channels, n_freqs).
#     n_bins: int
#         Number of neighbouring frequency bins to use to compute standard
#         deviation. Bins will be taken from each side at each frequency,
#         thus altogether 2*n_bins bins will be used.
#     mode: str ('z' | 'baseline' | 'snr')
#         Whether to do subtract baseline ('baseline'), divide by standard
#         deviation ('snr'), or both ('z', first baseline then snr).
#         Defaults to 'z'.
#     n_gap: int
#         Gap between target frequency and neighbouring bins.
#     minmax: Bool
#         Whether or not to remove minimum and maximum values from z-scoring.
#     Returns:
#     psds_snr: array
#         The PSD as SNRs, shape (n_channels, n_freqs).
#         SNR is amplitude divided by standard deviation.
#         If z-score required, subtract mean separately using baseline
#         correction.
#     """
#     # TO DO: faster with matrix computations?
#     psds_z = deepcopy(psds)  # initialise output array

#     for ff in np.arange(0, psds.shape[1]):  # for frequencies

#         # take care of edges in PSD
#         m = np.max([0, ff - n_bins - n_gap])
#         n = np.min([psds.shape[1], ff + 1 + n_bins + n_gap])

#         for cc in np.arange(0, psds.shape[0]):  # for channels

#             # neighbouring elements before and after this frequency
#             base_idx = np.r_[np.arange(m, ff - n_gap), np.arange(ff + 1 +
#                                                                  n_gap, n)]
#             # get baseline amplitudes
#             baseline = psds[cc, base_idx]

#             # indices of minimum and maximum baseline value
#             minidx = baseline.argmin()
#             maxidx = baseline.argmax()

#             # remove min/max values from baseline
#             np.delete(baseline, (minidx, maxidx))

#             if mode in ['z', 'baseline']:

#                 # baseline-correct at this frequency
#                 psds_z[cc, ff] = psds_z[cc, ff] - np.average(baseline)

#             if mode in ['z', 'snr']:

#                 # Compute SNR for one frequency
#                 psds_z[cc, ff] = psds_z[cc, ff] / np.std(baseline)

#     return psds_z


### now part of psd_z_score()
# def psd_correct_baseline(psds, n_bins, n_gap=0):
#     """BASELINE-CORRECT PSDs with neighbouring frequency bins.

#     Parameters:
#     psds: instance of Evoked
#         PSD, data of shape (n_channels x n_freqs),
#         from mne.time_frequency.psd_welch
#     n_gap: int
#         Gap between target frequency and neighbouring bins.

#     Returns:
#     psds_base: array
#         Baseline-corrected PSD, shape (n_channels x n_freqs)
#     """
#     data = psds.data

#     # initialise output
#     psds_base = deepcopy(data)

#     data_base = np.zeros(data.shape)  # initialise output array
#     for ff in np.arange(0, data.shape[1]):  # for frequencies

#         # take care of edges in PSD
#         m = np.max([0, ff - n_bins - n_gap])
#         n = np.min([data.shape[1], ff + n_bins + n_gap])

#         for cc in np.arange(0, data.shape[0]):  # for channels

#             # neighbouring elements before and after this frequency
#             base_idx = np.r_[np.arange(m, ff - n_gap), np.arange(ff + 1 +
#                                                                  n_gap, n)]

#             # get baseline amplitudes
#             baseline = data[cc, base_idx]

#             # baseline-correct at this frequency
#             data_base[cc, ff] = data[cc, ff] - np.average(baseline)

#     # insert data into Evoked object
#     psds_base.data = data_base

#     return psds_base


# Plot PSDs in Evoked format
def plot_psd_as_evo(psd, sbj_path, picks=None, txt_label='',
                    close_fig=1, scalings=dict(eeg=1e6, grad=1e13, mag=1e15)):
    """Plot PSD disguised as instance of Evoked.

    Parameters:
        psd: The PSD as Evoked object
        sbj_path: path where "Figures" sub-directory is
        picks: str | list | slice | None
            As for plot_joint:
            Channels to include. Slices and lists of integers will be
            interpreted as channel indices. In lists, channel type strings
            (e.g., ['meg', 'eeg']) will pick channels of those types, channel
            name strings (e.g., ['MEG0111', 'MEG2623'] will pick the given
            channels. Can also be the string values “all” to pick all channels,
            or “data” to pick data channels. None (default) will pick all
            channels.
        txt_label: string to make filename specific (at beginning)
        close_fig: close figure (1) or not (0)
        scalings: dict, scalings for eeg/grad/mag
    Returns:
        figs: list
        The list of pyplot figures created.
    """
    psd_as_evo = deepcopy(psd)  # will be modified

    # keep a copy for scaling below
    psd_tmp = deepcopy(psd_as_evo)

    # CROP PSD for display
    psd_as_evo.crop(tmin=config.crop_times[0], tmax=config.crop_times[1])

    # EEG present?
    is_eeg = psd_as_evo.__contains__('eeg')

    # quick hack to get scaling approximately right
    # avoid first 1Hz in scaling
    psd_tmp.crop(config.crop_times[0], tmax=config.crop_times[1])

    # Default mne-python scalings for evo.plot(), just for clarity
    # scalings = dict(eeg=1e6, grad=1e13, mag=1e15)
    # show y-axis with original values
    # scalings = dict(eeg=1., grad=1., mag=1.)

    # units for y-axis of PSD plots
    units = {'mag': r'amp/$\sqrt{Hz}$', 'grad': r'amp/$\sqrt{Hz}$'}

    ch_types = ['mag', 'grad']

    if is_eeg:  # only if EEG in fiff-file
        units['eeg'] = r'amp/$\sqrt{Hz}$'
        ch_types.append('eeg')

    # plot y-axis range (can be negative after baseline-correction)
    ylim = {'mag': [np.min(psd_tmp.data[2:306:3, :]) * scalings['mag'],
                    np.max(psd_tmp.data[2:306:3, :]) * scalings['mag']],
            'grad': [np.min([psd_tmp.data[0:306:3, :],
                            psd_tmp.data[1:306:3, :]]) * scalings['grad'],
                     np.max([psd_tmp.data[0:306:3, :],
                            psd_tmp.data[1:306:3, :]]) * scalings['grad']]}
    if is_eeg:
        ylim['eeg'] = [np.min(psd_tmp.data[306:376, :] * scalings['eeg']),
                       np.max(psd_tmp.data[306:376, :] * scalings['eeg'])]

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
                   ylim=ylim, time_unit='s')

    figs = []
    for ch_type in ch_types:

        # if channel selection specified, pick only channel names of the
        # current channel type
        if picks is not None:

            if ch_type == 'mag':

                picks_type = [chn for chn in picks if chn[-1] == '1']

            if ch_type == 'grad':

                picks_type = [chn for chn in picks if chn[-1] in ['2', '3']]

            if ch_type == 'eeg':

                picks_type = [chn for chn in picks if chn[:3] == 'EEG']

        else:

            picks_type = ch_type

        topomap_args = dict(scalings=scalings, time_format='%.2f Hz',
                            time_unit='ms', ch_type=ch_type)

        fig = psd_as_evo.plot_joint(times=ftimes, title=txt_label,
                                    ts_args=ts_args, picks=picks_type,
                                    topomap_args=topomap_args)

        figs.append(fig)

    for (cc, ff) in zip(ch_types, figs):

        for fig_format in ['.jpg']:  # jpg doesn't work, png does

            fname_fig = op.join(sbj_path, 'Figures', txt_label + '_' +
                                str(cc) + fig_format)
            print('Saving figure to %s' % fname_fig)

            # Save PSD figure
            ff.savefig(fname_fig)

        # close PSD figure if in QSUB
        if close_fig:
            plt.close(ff)

    return figs


def grand_average_evoked_arrays(evokeds):
    """Average data arrays across Evoked objects.

    Averages data arrays irrespective of channel names.

    Parameters:
    evokeds: list of instances of Evoked
        The data to average. The dimension of the data in each instance of
        Evoked must be the same. The data arrays will be averaged irrespective
        of channels names.

    Returns:
    gm_evoked: instance of Evoked
        The averaged data as instance of Evoked. The info from the first
        item in evokeds is used.
    """
    # dimension expected for all data arrays
    m, n = evokeds[0].data.shape

    info = evokeds[0].info

    datas = []  # will contain data arrays across list items

    for evo in evokeds:

        data = evo.data

        # check that dimensions are all the same
        if data.shape != (m, n):

            print('Dimensions for averaging do not match!')

            return

        datas.append(data)

    # average across data arrays
    gm_data = np.mean(datas, axis=0)

    tmin = evokeds[0].times[0]

    nave = len(evokeds)

    gm_evoked = EvokedArray(gm_data, info, tmin=tmin, nave=nave)

    return gm_evoked


def svd_per_channel_type(evokeds, idx=None):
    """Compute singular values from SVD per channel type.
    Parameters:
    evokeds: list of intances of Evoked
        The evoked data to use.
    idx: array-like
        The indices along the time axis of evoked to include in SVD.
        If none (default), use all samples.

    Returns:
    ss: list dict of numpy arrays
        The singular values for each evoked dataset per channel type.
    """

    ss = []

    # make sure there is a list of evokeds
    if type(evokeds) is not list:

        evokeds = [evokeds]

    for evoked in evokeds:

        if idx is None:

            idx = evoked.data.shape[1]

        ss.append({})  # will contain singular values for evokeds and types

        for cht in ['mag', 'grad', 'eeg']:

            if evoked.__contains__(cht):  # only do channel types in evoked

                evo = deepcopy(evoked)  # pick_channels will modify it

                # pick the current channel type
                if cht == 'eeg':

                    meg = False
                    eeg = True

                else:

                    meg = cht
                    eeg = False

                evo = evo.pick_types(meg=meg, eeg=eeg)

                # get data as numpy array for desired samples
                data = evo.data[:, idx]

                # SVD, singular values only
                s = np.linalg.svd(data, compute_uv=False)

                # keep singular values per channel type for this evoked
                ss[-1][cht] = s

    return ss


# NOT USED
# plot TRF results
def plot_TFR(powtfr, idx_std, std_max_channel, freq_std, sbj_path, raw_stem_in,
             txt_label, close_fig):
    """Time-frequency analysis of runs. Currently not used."""

    # powtfr: TFRAverage object
    # idx_std: index to standard frequency
    # std_max_channel: name of channel with max amplitude at standard frequency
    # freq_std: standard presentation frequency
    # sbj_path: path where "Figures" sub-directory is
    # raw_stem_in: part of the filename, for JPG filename
    # txt_label: string to make filename more specific
    # close_fig: close figure or not

    # get maximum value across time at standard frequency for scaling
    idx_fr = np.argmin(np.abs(config.tfr['freqs'] - freq_std))
    vmax = powtfr.data[idx_std, idx_fr, :].max()

    # plot TFR ratio for maximum channel determined above
    fig = powtfr.plot(picks=std_max_channel, title=std_max_channel, vmin=-vmax,
                      vmax=vmax)

    # yticks at integer frequencies, labels every 5Hz
    ticks = np.arange(int(config.tfr['freqs'][0]), int(config.tfr['freqs'][-1]), 1)
    lab_freqs = np.arange(int(config.tfr['freqs'][0]), int(config.tfr['freqs'][-1]), 5)
    labels = [x if x in lab_freqs else '' for x in ticks]
    plt.yticks(ticks=ticks, labels=labels)

    crop_str = str(int(10000. * config.crop_times[0])) + '_' + str(int(10000. * config.crop_times[1]))

    # name depends on crop_time in config.py
    fname_fig = op.join(sbj_path, 'Figures', 'TFR_' + txt_label + raw_stem_in + '_sss_f_raw_ica_' + crop_str + '.jpg')
    print('Saving figure to %s' % fname_fig)

    # Save PSD figure
    fig.savefig(fname_fig)

    # close PSD figure if in QSUB
    if close_fig:
        plt.close(fig)

#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
Get data segments for individual frequency sweeps for FPVS Frequency Sweep.

Read raw data, find runs, segment into individual frequency sweeps,
average sweeps across runs, write the average as raw file
(and ascii file if desired).
Compute TFR if specified.
==========================================

OH, October 2019
"""

import sys

from os import path as op
import numpy as np

from copy import deepcopy

# import matplotlib
# matplotlib.use('Agg') #  for running graphics on cluster ### EDIT

from matplotlib import pyplot as plt

from importlib import reload

# Code to save in EDF format
# https://gist.github.com/skjerns/bc660ef59dca0dbd53f00ed38c42f6be
from save_edf import write_edf

import mne

import config_sweep as config
reload(config)


print(mne.__version__)

# endings of raw files, e.g. with/without ICA
raw_fname_end = 'f_raw_ica.fif'
# raw_fname_end = 'f_raw.fif'

# perform TFR of raw data or not
do_tfr = config.do_tfr

close_fig = 1  # close figures only if close_fig==1

ascii_eeg_edf = 1  # write average raw data for EEG also as ascii and edf files?

# plt.ion() # interactive plotting


def run_get_sweeps(sbj_id):
    """Compute spectra for one subject."""

    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    fig = {}  # figures

    # initialise dict for results
    # one dict per condition (e.g. 'hflf'), then for frequency (e.g. '12.'),
    # then list of individual sweeps
    names = []  # names of conditions
    for raw_stem_in in sss_map_fname[1][2:]:

        names.append(raw_stem_in[:4])

    names = np.unique(names)

    data_runs = {}
    for name in names:

        data_runs[name] = {}

        if name == 'face':

            # faces only have one frequency
            data_runs[name]['6.0'] = []

        else:

            for freq in config.fpvs_freqs:

                data_runs[name][str(freq)] = []


    # Note: Skip first two files since they are rest, no events
    for raw_stem_in in sss_map_fname[1][2:]:

        # raw_fname_in = op.join(sbj_path, raw_stem_in[:-3] + 'f_raw_ica.fif')
        raw_fname_in = op.join(sbj_path, raw_stem_in[:-3] + raw_fname_end)

        print('\n###\nReading raw file %s.' % raw_fname_in)
        raw_ori = mne.io.read_raw_fif(raw_fname_in, preload=True)

        raw = deepcopy(raw_ori)  # keep raw_ori for possible TFR analysis

        # Check if EEG data present in fiff-file
        is_eeg = raw.__contains__('eeg')

        eeg, meg = is_eeg, True

        event_file = op.join(sbj_path, raw_stem_in + '_sss_f_raw-eve.fif')
        print('Reading events from %s.' % event_file)
        events = mne.read_events(event_file)

        # Find indices of good events (onsets of runs without missing frames)
        event_ids = config.fpvs_event_ids
        run_duration = config.fpvs_n_sweeps * config.fpvs_sweep_duration + \
            config.fpvs_leadin

        # idx_good, idx_bad: lists of indices to onsets of good/bad runs
        idx_good, idx_bad = find_good_events(events, event_ids=event_ids,
                                             run_duration=run_duration,
                                             sfreq=raw.info['sfreq'])

        print('Good runs:')
        print(idx_good)
        print(events[idx_good, :])

        if len(idx_bad) != 0:
            print('Bad runs:')
            print(idx_bad)
            print(events[idx_bad, :])
        else:
            print('No bad runs.')

        # go through all indices to good runs
        for idx in idx_good:

            # onset time (s) for this good run
            # there is one second gap between trigger and stimulus onset
            # note: samples don't start at 0, but times do
            onset_time = (events[idx, 0] - raw.first_samp) / raw.info['sfreq'] + 1.

            if raw_stem_in[:4] == 'face':  # faces don't have "sweeps"

                # just get one "sweep"
                raw_sweep = get_sweeps_from_raw(raw, t0=onset_time,
                                                sweep_duration=60.,
                                                n_sweeps=1)

                # for faces use whole run
                data_runs[raw_stem_in[:4]]['6.0'].append(raw_sweep[0])

            else:  # for frequency sweeps

                n_sweeps = len(config.fpvs_freqs)

                print('ID: %d, idx: %d, onset time: %f.' % (events[idx, 2],
                                                            idx, onset_time))

                # raw_sweeps: list of raw instances per data segments,
                # one per frequency
                raw_sweeps = get_sweeps_from_raw(raw, t0=onset_time,
                                                 sweep_duration=config.fpvs_sweep_duration,
                                                 n_sweeps=n_sweeps)

                for [fi, ff] in enumerate(config.fpvs_freqs):

                    data_runs[raw_stem_in[:4]][str(ff)].append(raw_sweeps[fi])

    # average raw files per condition and frequency across runs
    # write the result as raw fiff-file
    for cond in data_runs.keys():  # conditions

        for freq in data_runs[cond].keys():  # frequencies

            raw_avg = average_raws(data_runs[cond][freq])

            # remove dot from frequency string
            fname = 'rawavg_%s_%s.fif' % (cond, ''.join(freq.split('.')))
            fname_raw_out = op.join(sbj_path, fname)

            print('Writing average raw data to %s:' % fname_raw_out)

            raw_avg.save(fname_raw_out, overwrite=True)

            # if ascii and edf data for EEG requested as well
            if ascii_eeg_edf:

                # reduce to EEG only
                raw_avg_eeg = raw_avg.pick_types(meg=False, eeg=True)

                # ASCII
                fname = 'rawavg_%s_%s_eeg.asc' % (cond, ''.join(freq.split('.')))
                fname_asc_out = op.join(sbj_path, fname)

                data_out = raw_avg_eeg.get_data()

                print('Writing average ASCII data (%d, %d) to %s:'
                      % (data_out.shape[0], data_out.shape[1], fname_asc_out))

                np.savetxt(fname_asc_out, data_out)

                # EDF
                fname = 'rawavg_%s_%s_eeg.edf' % (cond, ''.join(freq.split('.')))
                fname_edf_out = op.join(sbj_path, fname)

                print('Writing average EDF data to %s:' % fname_edf_out)

                write_edf(raw_avg_eeg, fname_edf_out, overwrite=True)

    return data_runs


def find_good_events(events, event_ids, run_duration, sfreq):
    """Find the onsets of good runs in raw data.
    Parameters:
    events: events from raw data
    event_ids: list of int, possible triggers of run onsets
    run_duration: duration (s) of a run within session
    sfreq: float, sampling frequency (Hz)
    Returns:
    idx_good: list of indices to onsets of good runs
    idx_bad: list of indices to onsets of bad runs
    """
    max_missed = 2  # how many frames turn a run invalid

    idx_good, idx_bad = [], []  # initialise output

    # number of indices for events in this run
    n_idx = int(run_duration * sfreq)

    print(n_idx)

    # find all onsets in events based on event_ids
    onsets = [ee for ee in events if (ee[2] in event_ids)]

    for onset in onsets:

        # find index of this event
        onset_idx = np.where(events[:, 0] == onset[0])[0][0]

        # get all indices for events in this run
        idx_run = np.where((events[:, 0] > onset[0]) &
                           (events[:, 0] < onset[0] + n_idx))[0]

        print(idx_run[0])
        print(idx_run[-1])

        # get all events for this run
        events_run = events[idx_run, :]

        # check if missed frames present, and how many
        missed_frames = np.where(events_run[:, 2] == 20)[0]

        print(missed_frames)

        # if good run found
        if (len(missed_frames) == 0) or (missed_frames.shape[0] < max_missed):

            idx_good.append(onset_idx)

        else:  # if invalid due to missing frames

            idx_bad.append(onset_idx)

    return idx_good, idx_bad


def get_sweeps_from_raw(raw, t0, sweep_duration, n_sweeps):
    """Get segments from raw data for individual frequency sweeps.
    Parameters:
    raw: instance of Raw
    t0: start time of segment in s
    sweep_duration: duration of one sweep at one frequency (s)
    n_sweeps: int, number of sweeps (frequencies) per run
    Returns:
    raw_sweeps: list of raw instances per data segments, one per frequency
    """
    raw_sweeps = []  # initialise output

    for ss in np.arange(0, n_sweeps):

        tmin = t0 + (ss * sweep_duration)
        tmax = t0 + (ss + 1) * sweep_duration

        raw_cp = raw.copy()

        raw_cp.crop(tmin=tmin, tmax=tmax)

        raw_sweeps.append(raw_cp)

    return raw_sweeps


def average_raws(raws):
    """Average data across raw files.
    Parameters:
    raws: list of instances of Raw
    Returns:
    raw_avg: instance of Raw
    Contains average of raw data from raws.
    Every item of raws must have same info.
    """
    # get data array from first file
    data = raws[0].get_data()

    if len(raws) > 1:

        for raw in raws[1:]:

            data += raw.get_data()

        data = data / len(raws)

    # con't understand 'copy' option, using default
    raw_avg = mne.io.RawArray(data, raws[0].info, first_samp=0, copy='auto')

    return raw_avg





# get all input arguments except first
if len(sys.argv)==1:

    sbj_ids = np.arange(0,len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    # raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
    data_runs = run_get_sweeps(ss)

print('Done.')
#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Filter and clean data from FPVS Frequency Sweep.

Get event information, interpolate bad
EEG channels, (notch) filter.
==========================================

OH, modified October 2019
"""

import sys
from os import path as op
import numpy as np

# import matplotlib
# matplotlib.use('Agg')  # for running graphics on cluster ### EDIT
from matplotlib import pyplot as plt

from importlib import reload

import mne

import config_sweep as config
reload(config)

print('MNE Version: %s\n\n' % mne.__version__)  # just in case
print(mne)

# whether to show figures on screen or just write to file
show = False


def run_filter_raw(sbj_id):
    """Clean data for one subject."""
    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    bad_eeg = config.bad_channels[sbj_id]['eeg']  # bad EEG channels

    for raw_stem_in in sss_map_fname[1]:

        # input file to read
        raw_fname_in = op.join(sbj_path, raw_stem_in + '.fif')

        # result file to write
        raw_fname_out = raw_fname_in[:-7] + 'f_raw.fif'

        print('\n###\nReading raw file %s.' % raw_fname_in)

        raw = mne.io.read_raw_fif(raw_fname_in, preload=True)

        print('Fixing coil types.')
        raw.fix_mag_coil_types()

        # ONLY FOR EEG
        if any('EEG' in ch for ch in raw.info['ch_names']):

            print('Marking bad EEG channels: %s' % bad_eeg)
            raw.info['bads'] = bad_eeg

            print('Interpolating bad EEG channels.')
            raw.interpolate_bads(mode='accurate', reset_bads=True)

            print('Setting EEG reference.')
            raw.set_eeg_reference(ref_channels='average')

        else:

            print('No EEG channels found.\n')

        print('Applying Notch filter.')

        raw.notch_filter(np.array([50, 100]), fir_design='firwin')

        # str() because of None
        print('Applying band-pass filter %s to %s Hz.' % (str(config.l_freq),
                                                          str(config.h_freq)))

        # broad filter, including VGBR and ASSR frequencies
        # most settings are the MNE-Python defaults (zero-phase FIR)
        # https://mne.tools/dev/auto_tutorials/discussions/plot_background_filtering.html
        raw.filter(l_freq=config.l_freq, h_freq=config.h_freq, method='fir',
                   fir_design='firwin', filter_length='auto',
                   l_trans_bandwidth='auto', h_trans_bandwidth='auto')

        print('Saving data to %s.' % raw_fname_out)
        raw.save(raw_fname_out, overwrite=True)

        print('Finding events.')
        # note: short event duration
        events = mne.find_events(raw, stim_channel='STI101',
                                 consecutive='increasing', min_duration=0.002,
                                 verbose=True)

        # correct for stimulus presentation delay
        stim_delay = int(config.delay * raw.info['sfreq'])
        events[:, 0] = events[:, 0] + stim_delay

        event_file = op.join(sbj_path, raw_stem_in + '_sss_f_raw-eve.fif')
        print('Saving events to %s.' % event_file)
        mne.write_events(event_file, events)

        # plot only if events were found
        if events.size != 0:

            fig = mne.viz.plot_events(events, raw.info['sfreq'], show=show)

            fname_fig = op.join(sbj_path, 'Figures',
                                raw_stem_in + '_sss_f_raw_eve.pdf')
            print('Saving figure to %s' % fname_fig)

            fig.savefig(fname_fig)

            plt.close(fig)

        else:

            print('No events found in file %s.' % raw_fname_in)

    return raw, events


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    [raw, events] = run_filter_raw(ss)

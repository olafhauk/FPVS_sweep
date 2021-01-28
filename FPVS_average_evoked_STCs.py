#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Average source estimates of evoked data in FPVS.

Read morphed STCs produced in FPVS_source_estimation_evoked.py.
Average and plot them.
==========================================

OH, May 2020
"""

import sys

import os
from os import path as op
from os import stat

import numpy as np

# needed to run on SLURM
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from mayavi import mlab
mlab.options.offscreen = True

from copy import deepcopy

from importlib import reload

import mne
from mne.source_estimate import SourceEstimate

import config_sweep as config
reload(config)


print(mne.__version__)

# sub-directory for figures per subject
# separate for ICAed and non-ICAed data
if 'ica' in config.raw_ICA_suff:
    figs_dir = 'Figures_ICA'
else:
    figs_dir = 'Figures'

# output directory for figures
figs_path = op.join(config.grandmean_path, figs_dir)

# conditions
# conds = ['face', 'pwhf', 'pwlf', 'lfhf']
conds = config.do_conds

freqs_all = [str(ff) for ff in config.fpvs_freqs]


def run_average_STCs_evoked():
    sbj_ids = config.do_subjs

    # for evoked created with and without Notch filter for base frequency
    for do_notch in [0, 1]:

        if do_notch:  # if Notch filter at base frequency requested

            # add to epoch file name
            str_notch = '_nch'

        else:

            str_notch = ''

        for cond in conds:  # conditions

            if cond == 'face':  # hack, no frequency sweep for faces

                freqs = ['6.0']

            else:  # for all word condition, use all sweep frequencies

                freqs = freqs_all

            for freq in freqs:  # frequencies

                stcs = []

                for sbj_id in sbj_ids:

                    # path to subject's data
                    sbj_path = op.join(
                        config.data_path, config.map_subjects[sbj_id][0])

                    stc_fname = op.join(
                        sbj_path, 'STC', '%s_f_%s_%s%s_mph' %
                        (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                         str_notch))

                    print('Reading source estimate from %s.' % stc_fname)

                    stc = mne.read_source_estimate(stc_fname)

                    stcs.append(stc)

                print('Averaging %d STCs.' % len(stcs))
                avg_data = np.average([s.data for s in stcs], axis=0)

                print('0')

                # turn average into source estimate object
                stc_avg = SourceEstimate(avg_data, stcs[0].vertices,
                                         stcs[0].tmin, stcs[0].tstep)

                print('1')

                # path to subject's data
                sbj_path = op.join(config.data_path, 'GM')

                stc_fname = op.join(
                    sbj_path, 'STC', '%s_f_%s_%s%s' %
                    (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                     str_notch))

                times_plot = [0.13, 0.19, 0.24]

                print('2')

                for tt in times_plot:

                    # get index for time point to plot
                    idx = stc.time_as_index(times=[tt])

                    data = stc_avg.data[:, idx]

                    # threshold for plotting
                    thresh = np.abs(data).max()

                    # plot for left and right hemisphere
                    for hemi in ['both']:

                        # for some reason, 'both' only works for 'ven' but not
                        # for 'lat'
                        for view in ['lat', 'ven', 'cau']:

                            # continue until figure output has decent size,
                            # i.e. it didn't fail
                            st_size = 0

                            max_size = 10000  # bytes

                            max_tries = 3  # try 3 times

                            tries = 1

                            while st_size < max_size and tries <= max_tries:

                                if tries > 1:

                                    print('Trying again (%d).' % tries)

                                brain = stc_avg.plot(
                                    subject='fsaverage', initial_time=tt,
                                    time_label='',
                                    subjects_dir=config.subjects_dir,
                                    clim=dict(kind='value',
                                              lims=[0, thresh / 2., thresh]),
                                    hemi=hemi, views=view
                                )

                                # time in ms for file name
                                time_str = str(int(1000 * tt))

                                fname_fig = op.join(
                                    figs_path, '%s_f_%s_%s_%s_%s_%s%s.jpg' %
                                    (cond, config.raw_ICA_suff,
                                     ''.join(freq.split('.')), time_str,
                                     view, hemi, str_notch))

                                print('Saving figure to %s.' % fname_fig)

                                mlab.savefig(fname_fig)

                                st_size = stat(fname_fig).st_size

                                print('Size: %d' % st_size)

                                tries = tries + 1

                                mlab.close(all=True)

    return

# raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
run_average_STCs_evoked()

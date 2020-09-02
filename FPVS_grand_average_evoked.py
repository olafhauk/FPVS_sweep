#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Grand-average evoked data for FPVS.

Read evoked data created by FPVS_average_epochs.py,
average across subjects, plot curves and topographies.
==========================================

OH, May 2020
"""

### TO DO: maxfilter interpolate to common sensor configuration

import sys

import os
from os import path as op
import numpy as np

os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from mayavi import mlab
mlab.options.offscreen = True

import matplotlib
matplotlib.use('Agg')  #  for running graphics on cluster ### EDIT
from matplotlib import pyplot as plt

from copy import deepcopy

from importlib import reload

import mne

import config_sweep as config
reload(config)


print(mne.__version__)

# conditions
conds = config.do_conds

# base frequencies as strings
freqs_all = [str(ff) for ff in config.fpvs_freqs]

# sub-directory for figures per subject
# separate for ICAed and non-ICAed data
if 'ica' in config.raw_ICA_suff:
    figs_dir = 'Figures_ICA'
else:
    figs_dir = 'Figures'

# conditions
conds = config.do_conds

def run_grand_average_evoked(sbj_ids):
    """Plot evoked data for one subject."""
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

                evos = []

                for sbj_id in sbj_ids:

                    # path to subject's data
                    sbj_path = op.join(
                        config.data_path, config.map_subjects[sbj_id][0])

                    evo_fname = op.join(
                        sbj_path, 'AVE', '%s_f_%s_%s%s-ave.fif' %
                        (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                         str_notch))

                    print('Reading evoked data from %s.' % evo_fname)
                    # there is only one Evoked object in there
                    evoked = mne.read_evokeds(evo_fname, 0)

                    evos.append(evoked)

                # grand-average evoked data
                print('Grand-averaging %d files.' % len(evos))
                gm_evoked = mne.grand_average(evos)

                gm_fname = op.join(
                    config.grandmean_path, 'AVE', '%s_f_%s_%s%s-ave.fif' %
                    (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                     str_notch))

                print('Writing GM evoked to %s.' % gm_fname)

                mne.write_evokeds(gm_fname, gm_evoked)

                # parameters for plotting curves
                ts_args = dict(spatial_colors=True, gfp=True)

                # While we are here, plot evoked
                times = [0.12, 0.19, 0.24, 0.36]
                figs = gm_evoked.plot_joint(times=times, title=freq,
                                            ts_args=ts_args)

                # path to GM
                gm_path = op.join(
                    config.data_path, 'GM')

                # path to sub-directory for figures
                figs_path = op.join(gm_path, figs_dir)

                for [fi, fig] in enumerate(figs):

                    fig_fname = op.join(
                        figs_path, '%s_f_%s_%s%s_joint%d.jpg' %
                        (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                         str_notch, fi))

                    print('Saving figure to %s.' % fig_fname)

                    fig.savefig(fig_fname)

                plt.close('all')

    return


# get all input arguments except first
# if number not in config list, do all of them
if ((len(sys.argv) == 1) or
   (int(sys.argv[1]) > np.max(list(config.map_subjects.keys())))):

    # IDs don't start at 0
    sbj_ids = config.do_subjs

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]

# raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
data_runs = run_grand_average_evoked(sbj_ids)

#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Plot evoked data from FPVS sweeps for ERP analysis.

Read evoked created by FPVS_average_epochs.py,
plot curves and topographies.
==========================================

OH, May 2020
"""

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

# sub-directory for figures per subject
# separate for ICAed and non-ICAed data
if 'ica' in config.raw_ICA_suff:
    figs_dir = 'Figures_ICA'
else:
    figs_dir = 'Figures'

# conditions
conds = config.do_conds

def run_plot_evoked(sbj_id):
    """Plot evoked data for one subject."""
    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    # base frequencies as strings
    freqs_all = [str(ff) for ff in config.fpvs_freqs]

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

                evo_fname = op.join(
                    sbj_path, 'AVE', '%s_f_%s_%s%s-ave.fif' %
                    (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                     str_notch))

                print('Reading evoked data from %s.' % evo_fname)
                # there is only one Evoked object in there
                evoked = mne.read_evokeds(evo_fname, 0)

                # fig = evoked.plot(spatial_colors=True, gfp=True, time_unit='s')

                # fig_fname = evo_fname.replace('fif', 'jpg')

                # fig.savefig(fig_fname)

                # plotting parameters for plot_joint()

                # ts_args = dict(spatial_colors=True, scalings=scalings, units=units,
                #        ylim=ylim, time_unit='s')

                # topomap_args = dict(scalings=scalings, time_format='%.2f Hz',
                #                     time_unit='ms', ch_type=ch_type)

                figs = evoked.plot_joint(times='peaks', title=freq)

                # path to sub-directory for figures
                figs_path = op.join(sbj_path, figs_dir)

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
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    # raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
    data_runs = run_plot_evoked(ss)

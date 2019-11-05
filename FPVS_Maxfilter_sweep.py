#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
Maxfilter data from FPVS with Frequency Sweep.

=============================================
OH, modified October 2019
"""
import sys
import os
from os import path as op
import numpy as np

from importlib import reload

import config_sweep as config
reload(config)

MF = config.MF  # Maxfilter parameters


def run_maxfilter(sbj_id):
    """Run maxfilter for one subject."""
    # path to raw data for maxfilter
    map_subject = config.map_subjects[sbj_id]

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    raw_path = op.join(config.cbu_path, map_subject[0], map_subject[1])

    # use sixth filename for trans
    raw_fname_ref = op.join(raw_path,
                            sss_map_fname[0][config.MF['trans']] + '.fif')

    # maxfilter option for bad channels
    if config.bad_channels[sbj_id]['meg'] != []:

        # bad MEG channels without 'MEG'
        bad_channels = config.bad_channels[sbj_id]['meg']

        bads = [chn[3:] for chn in bad_channels]

        bad_cmd = '-bad %s' % ' '.join(bads)

    else:

        bad_cmd = ''

    for (raw_fname_in, raw_fname_out) in zip(sss_map_fname[0],
                                             sss_map_fname[1]):

        fname_in = op.join(config.cbu_path, map_subject[0], map_subject[1],
                           raw_fname_in + '.fif')

        fname_out = op.join(config.data_path, map_subject[0],
                            raw_fname_out + '.fif')

        log_fname_out = op.join(config.data_path, map_subject[0],
                                raw_fname_out + '_sss.log')  # log-file

        if op.exists(fname_out):  # if necessary delete existing MF output file
            os.remove(fname_out)

        # replace period to avoid confusion with file names
        if MF['st_duration'] is None:
            st_cmd = ''
        else:
            st_cmd = ' -st %s -corr %f' % (str(int(MF['st_duration'])),
                                           MF['st_correlation'])

        origin = MF['origin']
        ori_cmd = ' -origin %.0f %.0f %.0f ' % (1000 * origin[0], 1000 *
                                                origin[1], 1000 * origin[2])

        order_cmd = '-in %d  -out %d' % (MF['in'], MF['out'])

        if not(MF['mv'] == ''):
            mv_cmd = '-movecomp %s' % MF['mv']

        mf_cmd = '  %s \
                    -f %s \
                    -o %s \
                    -trans %s \
                    -frame %s \
                    -regularize %s \
                    %s \
                    %s \
                    %s \
                    %s \
                    %s \
                    -autobad on \
                    -force \
                    -linefreq 50 \
                    -v \
                    | tee %s' \
                    % (MF['NM_cmd'],
                       fname_in,
                       fname_out,
                       raw_fname_ref,
                       MF['frame'],
                       MF['regularize'],
                       st_cmd,
                       ori_cmd,
                       order_cmd,
                       mv_cmd,
                       bad_cmd,
                       log_fname_out)

        print('Maxfilter command: %s' % mf_cmd)

        # Execute maxfilter command
        os.system(mf_cmd)

# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    run_maxfilter(ss)

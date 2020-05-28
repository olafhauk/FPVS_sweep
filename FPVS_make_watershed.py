#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
=========================================================
Create BEM surfaces for FPVS using watershed algorithm.

Needs Freesurfer. Run freesurfer_6.0.0 before starting
mne_python_v0.19.
=========================================================

"""
# OH, March 2020

print(__doc__)

import sys

import os
from os import path as op

import numpy as np

import importlib
from importlib import reload

import mne

import config_sweep as config
reload(config)

# # For MNE and Freesurfer
# os.environ = {
# 'MNE_ROOT': '/imaging/local/software/mne/mne_2.7.3/x86_64/MNE-2.7.3-3268-Linux-x86_64/',
# 'MNI_PERL5LIB': '/imaging/local/software/freesurfer/6.0.0/x86_64/mni/share/perl5',
# 'MNI_DATAPATH': '/imaging/local/software/freesurfer/6.0.0/x86_64/mni/data',
# 'FMRI_ANALYSIS_DIR': '/imaging/local/software/freesurfer/6.0.0/x86_64/fsfast',
# 'FREESURFER_HOME': '/imaging/local/software/freesurfer/6.0.0/x86_64',
# 'LD_LIBRARY_PATH': '/imaging/local/linux/lib:/imaging/local/software/freesurfer/6.0.0/x86_64/lib_flat',
# 'PERL5LIB': '/imaging/local/software/freesurfer/6.0.0/x86_64/mni/share/perl5:/imaging/local/lib/perl5',
# 'MNI_DIR': '/imaging/local/software/freesurfer/6.0.0/x86_64/mni',
# 'LOCAL_DIR': '/imaging/local/software/freesurfer/6.0.0/x86_64/local',
# 'PATH': '/imaging/local/software/freesurfer/6.0.0/x86_64/bin:/imaging/local/software/freesurfer/6.0.0/x86_64/fsfast/bin:/imaging/local/software/freesurfer/6.0.0/x86_64/tktools:/imaging/local/software/freesurfer/6.0.0/x86_64/mni/bin:/hpc-software/bin:/imaging/local/software:/opt/TurboVNC/bin:/opt/TurboVNC/bin:/imaging/local/software/anaconda/latest/x86_64/bin:/opt/gold/bin:/opt/gold/sbin:/hpc-software/bin:/imaging/local/software:/cluster-software/gold/2.2.0.5/sbin:/cluster-software/maui/bin:/cluster-software/torque-2.3/sbin:/cluster-software/torque-2.3/bin:/opt/TurboVNC/bin:/opt/TurboVNC/bin:/usr/lib64/qt-3.3/bin:/imaging/local/software/anaconda/latest/x86_64/bin:/opt/gold/bin:/opt/gold/sbin:/usr/local/bin:/bin:/usr/bin:/opt/bin:/usr/local/wrapper/bin:/usr/local/bin:/usr/bin/X11:/usr/local/cbu/scripts/:/imaging/local/linux/bin:/imaging/local/linux/bin/elekta:/imaging/local/linux/bin/tractography/bin:/imaging/local/scripts:/home/olaf/bin:/cbu/language/celex/bin:./:/usr/local/wrapper/bin:/usr/local/bin:/usr/bin/X11:/usr/local/cbu/scripts/:/imaging/local/linux/bin:/imaging/local/linux/bin/elekta:/imaging/local/linux/bin/tractography/bin:/imaging/local/scripts:/home/olaf/bin:/cbu/language/celex/bin:./',
# 'MINC_BIN_DIR': '/imaging/local/software/freesurfer/6.0.0/x86_64/mni/bin',
# 'FSFAST_HOME': '/imaging/local/software/freesurfer/6.0.0/x86_64/fsfast',
# 'MINC_LIB_DIR': '/imaging/local/software/freesurfer/6.0.0/x86_64/mni/lib',
# }

# Failed attempt to send output to both a file and stderr
import logging
logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())

# Freesurfer and MNE environment variables
filename = "/imaging/local/software/mne_python/set_MNE_2.7.3_FS_6.0.0_environ.py"
# for Python 3 instead of execfile
exec(compile(open(filename, "rb").read(), filename, 'exec'))


# where MRIs are
os.environ['SUBJECTS_DIR'] = config.subjects_dir


def run_make_watershed(sbj_id):

    subject = config.mri_subjects[sbj_id]

    if subject == '':

        print('No subject name for MRI specified - doing nothing now.')

        return

    print('Making BEM surfaces for %s.' % subject)

    print('Creating BEM surfaces.')

    # print(os.environ.keys())

    # mne_cmd = '/imaging/local/software/mne/mne_2.7.3/x86_64/MNE-2.7.3-3268-Linux-x86_64//bin/mne_watershed_bem \
    #           --subject %s --overwrite' % subject

    # os.system(mne_cmd)

    mne.bem.make_watershed_bem(subject, subjects_dir=config.subjects_dir, atlas=True,
                               overwrite=True, copy=True)

# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]

for ss in sbj_ids:

    run_make_watershed(ss)

print('Done.')

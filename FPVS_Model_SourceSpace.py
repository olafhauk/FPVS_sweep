"""
=========================================================
BEM Model and Source Space for FPVS.
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print __doc__

# Russell's addition
import sys
# # for qsub: check matplotlib.use('Agg'), plt.ion(), plt.show(), do_show
# sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/sklearn/')
# sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/pysurfer/')
# sys.path.insert(1, '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/nibabel/')
# sys.path.insert(1, '/imaging/local/software/mne_python/v0.11/')
import matplotlib
matplotlib.use('Agg') # possibly for running on cluster

from importlib import reload

#sys.path.append('/imaging/local/software/python_packages/nibabel/2.0.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
# End

# Failed attempt to send output to both a file and stderr
import logging
logger = logging.getLogger()
logger.addHandler(logging.StreamHandler())

import mne
# mne.set_log_file(fname='/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/BEM_Model_and_SourceSpace.log', overwrite=None)

# where MRIs are
subjects_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/'

# where figures will be written to
bem_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/BEM_figs'

# subject IDs
subs = [
# me
'CBU160881',
# # her
# 'CBU000000',
# # the other one
# 'CBU000001',
]


if len(sys.argv)>1: # if in parallel mode
    print "Running subject(s) {0} now in PARALLEL mode".format(sys.argv)
    ss_idx = map(int, sys.argv[1:])
    subs_new = []
    for ii,ss in enumerate(ss_idx): # a bit cumbersome because lists cannot be used as indices
        subs_new.append(subs[ss])
    subs = subs_new
else:
    print "Running now in SERIAL mode"


conductivity_1 = (0.3,)  # for single layer
conductivity_3 = (0.3, 0.006, 0.3)  # for three layers

for ss in subs:
    print 
    print "###\nMaking BEM model for " + ss + "\n###"

    ### one-shell BEM for MEG
    print("MEG")
    model = mne.make_bem_model(subject=ss, ico=4,
                                conductivity=conductivity_1,
                                subjects_dir=subjects_dir)
    bem = mne.make_bem_solution(model)

    bem_fname = bem_dir + '/' + ss + '_MEG-bem.fif'

    print "###\nWriting BEM solution to " + bem_fname + "\n###"
    mne.bem.write_bem_solution(bem_fname, bem)

    ### three-shell BEM for EEG+MEG
    print('EEG+MEG')
    model = mne.make_bem_model(subject=ss, ico=4,
                                conductivity=conductivity_3,
                                subjects_dir=subjects_dir)
    bem = mne.make_bem_solution(model)

    bem_fname = bem_dir + '/' + ss + '_EEGMEG-bem.fif'

    print "###\nWriting BEM solution to " + bem_fname + "\n###"
    mne.bem.write_bem_solution(bem_fname, bem)


    ### set up source space
    print('Setting up source space.')
    src = mne.setup_source_space(ss, spacing='oct6',
                             subjects_dir=subjects_dir,
                             add_dist=False)

    src_fname = bem_dir + '/' + ss + '-src.fif'

    print "###\nWriting source spaces to " + src_fname + "\n###"
    mne.write_source_spaces(src_fname, src, overwrite=True)
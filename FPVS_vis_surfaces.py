"""
=========================================================
Visualise BEM surface for FPVS.
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print __doc__

# Russell's addition
import sys
# sys.path.append('/imaging/local/software/python_packages/nibabel/2.0.0')
# sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
# End

from importlib import reload

import mne


# where MRIs are
subjects_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/'

# where figures will be written to
fig_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/BEM_figs'

# subject IDs
subs = [
# me
'CBU160881',
# her
'CBU000000',
# the other one
'CBU000001',
]

### Function that plots BEM surfaces for one orientation
def plot_one_bem(ss, orientation, fig_dir):
    fig_bem = mne.viz.plot_bem(subject=ss, subjects_dir=subjects_dir, orientation=orientation, show=False)
    
    fig_fname = fig_dir + '/' + ss + '_' + orientation + '.pdf'
    print "Writing " + fig_fname
    fig_bem.savefig(fig_fname)

    return fig_bem

### Plot BEM surfaces one by one
orientations = ['coronal', 'axial', 'sagittal']

for ss in subs:
    print ss
    # plot BEM surfaces on MRI

    [plot_one_bem(ss, orientation, fig_dir) for orientation in orientations]
    




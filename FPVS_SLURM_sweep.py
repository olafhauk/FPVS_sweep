#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
==========================================
Submit sbatch jobs for FPVS Frequency Sweep
analysis
SLURM, Python 3
==========================================

OH, modified October 2019
"""

import subprocess
from os import path as op

print(__doc__)

# wrapper to run python script via qsub. Python3
fname_wrap = op.join('/', 'home', 'olaf', 'MEG', 'FPVS', 'MNE-Python',
                     'Python2SLURM.sh')

# indices of subjects to process
# subjs = range(1,5)
# EDIT
# subjs = [1, 2]
# subjs = [3, 4]
subjs = [1, 2, 3, 4]

job_list = [
    # # Neuromag Maxfilter
    # {'N':   'F_MF',                  # job name
    #  'Py':  'FPVS_Maxfilter_sweep',  # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '16G',                   # memory for qsub process
    #  'dep': '',                       # name of preceeding process (optional)
    #  'node': '--constraint=maxfilter'},  # node constraint for MF, just picked one
    # fix EEG electrode positions in fiff-files
    # # NOTE: Can get "Permission denied"
    # {'N':   'F_FE',                    # job name
    #  'Py':  'FPVS_fix_electrodes_sweep',      # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '1G',                    # memory for qsub process
    #  'dep': ''},                       # name of preceeding process (optional)
    # ### Filter raw data
    # {'N':   'F_FR',                  # job name
    #  'Py':  'FPVS_filter_raw_sweep',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '16G',                    # memory for qsub process
    #  'dep': ''},                      # name of preceeding process (optional)

    #  ### Compute ICA
    # {'N':   'F_CICA',                  # job name
    #  'Py':  'FPVS_Compute_ICA_sweep',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '4G',                    # memory for qsub process
    #  'dep': ''},                      # name of preceeding process (optional)
    #   ### Apply ICA
    # {'N':   'FP_AICA',                  # job name
    #  'Py':  'FPVS_Apply_ICA_sweep',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '2G',                    # memory for qsub process
    #  'dep': 'F_CICA'},                      # name of preceeding process (optional)


    ### Get sweeps from raw data and average
    {'N':   'F_GS',                  # job name
     'Py':  'FPVS_get_sweeps',          # Python script
     'Ss':  subjs,                    # subject indices
     'mem': '8G',                    # memory for qsub process
     'dep': ''},
    ### Compute PSDs for averaged sweeps and plot
    {'N':   'F_PSD',                  # job name
     'Py':  'FPVS_PSD_sweep',          # Python script
     'Ss':  subjs,                    # subject indices
     'mem': '8G',                    # memory for qsub process
     'dep': 'F_GS'},

    # ### Rename triggers
    # {'N':   'SR_RT',                  # job name
    #  'Py':  'SR_rename_triggers',     # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '4G',                    # memory for qsub process
    #  'dep': 'SR_FR'},                      # name of preceeding process (optional)
    # ### Perform ICA
    # {'N':   'SR_ICA',                 # job name
    #  'Py':  'SR_ICA_EOG',             # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '12G',                    # memory for qsub process
    #  'dep': 'SR_RT'},                      # name of preceeding process (optional)
    # ## Apply ICA
    # {'N':   'SR_AICA',                 # job name
    #  'Py':  'SR_apply_ICA_EOG',             # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '6G',                    # memory for qsub process
    #  'dep': 'SR_ICA'},                      # name of preceeding process (optional)
    # ### Epoch and average data
    # {'N':   'SR_EPO',                 # job name
    #  'Py':  'SR_Epoch',               # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '16G',                    # memory for qsub process
    #  'dep': ''},                      # name of preceeding process (optional)

    # # CHECK
    # # ### Epoch and average data - grating onset
    # # {'N':   'SR_EPOS',                 # job name
    # #  'Py':  'SR_Epoch_stim',               # Python script
    # #  'Ss':  subjs,                    # subject indices
    # #  'mem': '8G',                    # memory for qsub process
    # #  'dep': 'SR_EPO'},                      # name of preceeding process (optional)

    # ### compute TFR using wavelets
    # {'N':   'SR_TFW',                 # job name
    #  'Py':  'SR_TFR_Wavelet',         # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '10G',                    # memory for qsub process
    #  'dep': ''},                      # name of preceeding process (optional)
    # # ### compute TFR using wavelets - grating onset
    # {'N':   'SR_TFWS',                 # job name
    #  'Py':  'SR_TFR_Wavelet_stim',         # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '10G',                    # memory for qsub process
    #  'dep': 'SR_EPOS'},                      # name of preceeding process (optional)
    # ### Plot results from SR_TFR_Wavelet
    # {'N':   'SR_PTFW',                  # job name
    #  'Py':  'SR_Plot_TFR_Wavelet',      # Python script
    #  'Ss':  subjs,                      # subject indices
    #  'mem': '4G',                      # memory for qsub process
    #  'dep': ''},                        # name of preceeding process (optional)
    # #  ### Plot results from SR_TFR_Wavelet - grating onset
    # {'N':   'SR_PTFWS',                  # job name
    #  'Py':  'SR_Plot_TFR_Wavelet_stim',      # Python script
    #  'Ss':  subjs,                      # subject indices
    #  'mem': '4G',                      # memory for qsub process
    #  'dep': 'SR_TFWS'}                        # name of preceeding process (optional)
]

### Other processing steps
# ### compute band-limited time courses using Hilbert transform
    # {'N':   'SR_TFH',                 # job name
    #  'Py':  'SR_TFR_Hilbert',         # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '4G',                    # memory for qsub process
    #  'dep': 'SR_FR'},                      # name of preceeding process (optional)


# directory where python scripts are
dir_py = op.join('/', 'home', 'olaf', 'MEG', 'FPVS', 'MNE-Python')

# directory for qsub output
dir_sbatch = op.join('/', 'home', 'olaf', 'MEG', 'FPVS', 'MNE-Python',
                     'sbatch_out')


# keep track of qsub Job IDs
Job_IDs = {}

for job in job_list:

    for Ss in job['Ss']:

        Ss = str(Ss)  # turn into string for filenames etc.

        N = Ss + job['N'] # add number to front
        Py = op.join(dir_py, job['Py'])
        Cf = ''  # config file not necessary for FPVS
        mem = job['mem']

        # files for qsub output
        file_out = op.join(dir_sbatch,
                           job['N'] + '_' + Cf + '-%s.out' % str(Ss))
        file_err = op.join(dir_sbatch,
                           job['N'] + '_' + Cf + '-%s.err' % str(Ss))

        # if job dependent of previous job, get Job ID and produce command
        if 'dep' in job:  # if dependency on previous job specified
            if job['dep'] == '':
                dep_str = ''
            else:
                job_id = Job_IDs[Ss + job['dep'], Ss]
                dep_str = '--dependency=afterok:%s' % (job_id)
        else:
            dep_str = ''

        if 'node' in job:  # if node constraint present (e.g. Maxfilter)
            node_str = job['node']
        else:
            node_str = ''

        if 'var' in job:  # if variables for python script specified
            var_str = job['var']
        else:
            var_str = ''

        # sbatch command string to be executed
        sbatch_cmd = 'sbatch \
                        -o %s \
                        -e %s \
                        --export=pycmd="%s.py %s",subj_idx=%s,var=%s \
                        --mem=%s -t 1-00:00:00 %s -J %s %s %s' \
                        % (file_out, file_err, Py, Cf, Ss, var_str, mem,
                           node_str, N, dep_str, fname_wrap)

        # format string for display
        print_str = sbatch_cmd.replace(' ' * 25, '  ')
        print('\n%s\n' % print_str)

        # execute qsub command
        proc = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, shell=True)

        # get linux output
        (out, err) = proc.communicate()

        # keep track of Job IDs from sbatch, for dependencies
        Job_IDs[N, Ss] = str(int(out.split()[-1]))

#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
==========================================
Submit sbatch jobs for FPVS Frequency Sweep
analysis
SLURM, Python 3
==========================================

OH, modified October 2019
modified by Federica M for more subjects (ERP drive, MEG/FPVS/Scripts_Federica),
then re-adapted by OH Jan 2020
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
# subjs = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# subjs = [10, 11, 12, 13, 14] # up to FR excpet 13
# subjs = [15] 
# subjs = [16] up to Maxfilter
# subjs = [17]
# subjs = [13] 
# subjs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14] # done GS and PSD -> missing 13
# subjs = [13, 16, 17]
# subjs = [1] # check how looks after changing channels
# subjs = [1, 13, 15]
subjs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 15, 17, 18]
# subjs = [8]


job_list = [
    # Neuromag Maxfilter
    {'N':   'F_MF',                  # job name
     'Py':  'FPVS_Maxfilter_sweep',  # Python script
     'Ss':  subjs,                    # subject indices
     'mem': '16G',                   # memory for qsub process
     'dep': '',                       # name of preceeding process (optional)
     'node': '--constraint=maxfilter'},  # node constraint for MF, just picked one

    # # fix EEG electrode positions in fiff-files
    # # NOTE: Can get "Permission denied"; should be run separately
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

    # # ### Filter raw data FM-> generating also txt event file
    # # {'N':   'F_FR',                  # job name
    # #  'Py':  'FPVS_filter_raw_sweep_txtfile',          # Python script
    # #  'Ss':  subjs,                    # subject indices
    # #  'mem': '16G',                    # memory for qsub process
    # #  'dep': 'F_FR'},                      # name of preceeding process (optional)

    # #  ### Compute ICA
    # {'N':   'F_CICA',                  # job name
    #  'Py':  'FPVS_Compute_ICA_sweep',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '32G',                    # memory for qsub process
    #  'dep': 'F_FR'},                      # name of preceeding process (optional)
    #   ### Apply ICA (change ica_suff in config_sweep.py if necessary)
    # {'N':   'F_AICA',                  # job name
    #  'Py':  'FPVS_Apply_ICA_sweep',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '2G',                    # memory for qsub process
    #  'dep': 'F_CICA'},                      # name of preceeding process (optional)


    # ## Get sweeps from raw data and average (change ica_suff in config_sweep.py if necessary)
    # {'N':   'F_GS',                  # job name
    #  'Py':  'FPVS_get_sweeps',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '8G',                    # memory for qsub process
    #  'dep': 'F_AICA'},

    # ### Compute PSDs for averaged sweeps and plot (change ica_suff in config_sweep.py if necessary)
    # {'N':   'F_P_C',                  # job name
    #  'Py':  'FPVS_PSD_sweep_compute',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '8G',                    # memory for qsub process
    #  'dep': 'F_GS'},
    # ### Plot PSD results
    # {'N':   'F_P_P',                  # job name
    #  'Py':  'FPVS_PSD_sweep_plot',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '8G',                    # memory for qsub process
    #  'dep': 'F_P_C'},

    # ### Create Source Spaces
    # {'N':   'FP_SP',                  # job name
    #  'Py':  'FPVS_make_SourceSpace',          # Python script
    #  'Ss':  subjs,                    # subject indices
    #  'mem': '2G',                    # memory for qsub process
    #  'dep': ''}                      # name of preceeding process (optional)
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

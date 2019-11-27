#!/bin/bash
#PBS -q compute
#PBS -l walltime=24:00:00,mem=64GB
#PBS -o /home/olaf/MEG/FPVS/MNE-Python/qsub_out/reconall.out 
#PBS -e /home/olaf/MEG/FPVS/MNE-Python/qsub_out/reconall.err

# submit e.g. as qsub -t 0-21 FPVS_MRI_process.sh
# note: indices start at 0

# Runs mri_convert and recon-all to preprocess structural MRI data

export FSVER='5.3.0'

export FSDIR=${FSROOT}/${FSVER}

export FREESURFER_HOME=/imaging/local/software/freesurfer/${FSVER}/`arch`
export SUBJECTS_DIR=/imaging/`whoami`/subjects
export FUNCTIONALS_DIR=/imaging/`whoami`/sessions
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FREESURFER_HOME}/lib_flat

echo $FREESURFER_HOME
echo $SUBJECTS_DIR

source $FREESURFER_HOME/FreeSurferEnv.sh


export SUBJECTS_DIR=/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/
echo $SUBJECTS_DIR

subs=(\
# me
'CBU160881' \
# her
'CBU000000' \
# the other one
'CBU000001' \
)

if [ $PBS_ARRAYID == 0 ] # me
then
    echo "Now processing" ${subs[$PBS_ARRAYID]}

    mkdir /group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/${subs[$PBS_ARRAYID]}/mri/orig -p

    mri_convert /mridata/cbu/${subs[$PBS_ARRAYID]}_MR16002/20160923_171616/Series_002_xpace_MPRAGE_32chn_wt2500_OFF/1.3.12.2.1107.5.2.43.67035.2016092317315218252708598.dcm \
     /group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/${subs[$PBS_ARRAYID]}/mri/orig/001.mgz

    recon-all -subjid ${subs[$PBS_ARRAYID]} -all

elif [ $PBS_ARRAYID == 1 ] # her
then
    echo "Now processing" ${subs[$PBS_ARRAYID]}

    mkdir /group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/${subs[$PBS_ARRAYID]}/mri/orig -p

    mri_convert /imaging/gr02/ForOlaf/002/T1.nii.gz \
     /group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/${subs[$PBS_ARRAYID]}/mri/orig/001.mgz

    recon-all -subjid ${subs[$PBS_ARRAYID]} -all

elif [ $PBS_ARRAYID == 2 ] # the other one
then
    echo "Now processing" ${subs[$PBS_ARRAYID]}

    mkdir /group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/${subs[$PBS_ARRAYID]}/mri/orig -p

    mri_convert /imaging/gr02/ForOlaf/001/T1.nii \
     /group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/${subs[$PBS_ARRAYID]}/mri/orig/001.mgz

    recon-all -subjid ${subs[$PBS_ARRAYID]} -all

else
  echo "Unrecognised subject number"
fi 
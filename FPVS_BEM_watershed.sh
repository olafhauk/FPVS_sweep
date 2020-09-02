#!/bin/bash
#PBS -q compute
#PBS -l walltime=12:00:00,mem=64GB
#PBS -o /home/olaf/MEG/FPVS/MNE-Python/qsub_out/watershed.out 
#PBS -e /home/olaf/MEG/FPVS/MNE-Python/qsub_out/watershed.err

# submit e.g. as qsub -t 0-21 FPVS_BEM_watershed.sh -o ./qsub_out/BEM_watershed.out -e ./qsub_out/BEM_watershed.err 
# note: indices start at 0

# NOTE: watershed command may be commented out below

export FSVER='5.3.0'

export FSDIR=${FSROOT}/${FSVER}

export FREESURFER_HOME=/imaging/local/software/freesurfer/${FSVER}/`arch`
export SUBJECTS_DIR=/imaging/`whoami`/subjects
export FUNCTIONALS_DIR=/imaging/`whoami`/sessions
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FREESURFER_HOME}/lib_flat

echo $FREESURFER_HOME

source $FREESURFER_HOME/FreeSurferEnv.sh

export MNE_ROOT=/imaging/local/software/mne/mne_2.7.3/x86_64/MNE-2.7.3-3268-Linux-x86_64
export MNE_BIN_PATH=$MNE_ROOT/bin

export PATH=${PATH}:${MNE_BIN_PATH}
# source $MNE_ROOT/bin/mne_setup

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


export SUBJECT=${subs[$PBS_ARRAYID]}
echo $SUBJECT


# Computes the surfaces
$MNE_ROOT/bin/mne_watershed_bem --overwrite

# creates links in directory for MNE-Python:
# NOTE: unlink
ln -sf $SUBJECTS_DIR/$SUBJECT/bem/watershed/$SUBJECT'_inner_skull_surface' $SUBJECTS_DIR/$SUBJECT/bem/inner_skull.surf
ln -sf $SUBJECTS_DIR/$SUBJECT/bem/watershed/$SUBJECT'_outer_skull_surface' $SUBJECTS_DIR/$SUBJECT/bem/outer_skull.surf
ln -sf $SUBJECTS_DIR/$SUBJECT/bem/watershed/$SUBJECT'_outer_skin_surface'  $SUBJECTS_DIR/$SUBJECT/bem/outer_skin.surf
ln -sf $SUBJECTS_DIR/$SUBJECT/bem/watershed/$SUBJECT'_brain_surface'       $SUBJECTS_DIR/$SUBJECT/bem/brain_surface.surf

echo $SUBJECT " done"
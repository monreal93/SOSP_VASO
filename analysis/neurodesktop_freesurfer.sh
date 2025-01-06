cd /neurodesktop-storage/5T3/Alejandro/sosp_vaso/data


folder="08302024_sb_9T"
subject="t1_test1"
# scan="s2_brain_t1"
scan="sb_22_OUT_2shot_6te_14fa"
recon_all=0


cd ${folder}
cd ./analysis/
mkdir ./freesurfer

export SINGULARITYENV_SUBJECTS_DIR=./freesurfer

# Neurodesktop tools to use:
ml freesurfer
ml afni

# # To make "orig.mgz" match the proper orientation, I need to do this:
# mri_convert ./freesurfer/cv_01/mri/orig.mgz ./tmp/sample.nii
# 3dresample -orient RPS -prefix ./tmp/sample.nii -input ./tmp/sample.nii -overwrite
# 3drefit -orient LPI./tmp/sample.nii -overwrite
# 3dLRflip -Y -prefix ./tmp/sample.nii ./tmp/sample.nii

# t1_file=./${scan}_t1_full_brain.nii
t1_file=./${scan}/${scan}_t1.nii

if [ $recon_all -eq 1 ]; then
    #### ) recon-all
    recon-all -subject ${subject} -i ${t1_file} -all -notal-check
else
    #### ) autorecon1
    recon-all -autorecon1 -notal-check -notalairach -subject ${subject} -i ${t1_file}

    # If skull srip did a bad job, add some control points.. and then run -normalization
    # mri_normalize -f ./freesurfer/cv_01/tmp/control.dat ./freesurfer/cv_01/mri/orig.mgz ./freesurfer/cv_01/mri/orig.mgz
    # recon-all -normalization -subjid ${subject}

    # To Fix a bad skull stip, wsthresh values 1:50.. 
    # recon-all -skullstrip -wsthresh 5 -clean-bm -no-wsgcaatlas -subjid ${subject}

    # ### ) autorecon2
    # recon-all -autorecon2 -subject ${subject}

    # ### ) autorecon3
    # recon-all -autorecon3 -subject ${subject}
fi

# If WM segmentation fails, add some control points and thren recon again...
# https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/ControlPoints_freeview/
# recon-all -autorecon2-cp -autorecon3 -subjid cp_before

mkdir ./t1
# # saving WM and GM as masks..
mri_extract_label ./freesurfer/t1/mri/aseg.mgz 2 41 ./t1/wm_msk.nii
mri_extract_label ./freesurfer/t1/mri/aseg.mgz 3 42 ./t1/gm_msk.nii

mri_convert ./freesurfer/t1/mri/brain.mgz ./t1/brain.nii

# # Converting surface to gii or nii
# mris_convert ./freesurfer/ale_sv_orig/surf/lh.pial   ./tmp/lh_pial.gii
# mris_convert ./freesurfer/ale_sv_orig/surf/rh.pial   ./tmp/rh_pial.gii
# mris_convert ./freesurfer/ale_sv_orig/surf/lh.white  ./tmp/lh_white.gii
# mris_convert ./freesurfer/ale_sv_orig/surf/rh.white  ./tmp/rh_white.gii
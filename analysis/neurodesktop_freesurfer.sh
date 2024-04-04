cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data


folder="09192023_sv_abc_paper"
subject="t1"
# scan="registered_Warped_cv_01"
scan="my_brain_mp2rage_iso0.7mm_iPAT3_20230413161444_13"
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
# # 3dLRflip -Z -prefix ./tmp/sample.nii ./tmp/sample.nii -overwrite  # Not sure if I need this
# mri_convert ./tmp/orig.nii ./freesurfer/cv_01/mri/orig.mgz

# t1_file=./ants/${scan}.nii
t1_file=../raw/nifti/${scan}.nii

if [ $recon_all -eq 1 ]; then
    #### ) recon-all
    recon-all -subject ${subject} -i ${t1_file} -all -notal-check
else
    #### ) autorecon1
    recon-all -autorecon1 -notal-check -subject ${subject} -i ${t1_file}

    # If skull srip did a bad job, add some control points.. and then run -normalization
    # mri_normalize -f ./freesurfer/cv_01/tmp/control.dat ./freesurfer/cv_01/mri/orig.mgz ./freesurfer/cv_01/mri/orig.mgz
    # recon-all -normalization -subjid ${subject}

    # To Fix a bad skull stip, wsthresh values 1:50.. 
    # recon-all -skullstrip -wsthresh 35 -clean-bm -no-wsgcaatlas -subjid ${subject}

    #### ) autorecon2
    # recon-all -autorecon2 -subject ${subject}

    #### ) autorecon3
    # recon-all -autorecon3 -subject ${subject}
fi

# If WM segmentation fails, add some control points and thren recon again...
# https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/ControlPoints_freeview/
# recon-all -autorecon2-cp -autorecon3 -subjid cp_before

# # saving WM and GM as masks..
# mri_extract_label ./freesurfer/t1/mri/aseg.mgz 2 41 ./t1/wm_mask.nii
# mri_extract_label ./freesurfer/t1/mri/aseg.mgz 3 42 ./t1/gm_mask.nii

# # Converting surface to gii or nii
# mris_convert ./freesurfer/cv_01_full_brain/surf/lh.pial.T1 ./tmp/lh_pial.gii
# mris_convert ./freesurfer/cv_01_full_brain/surf/rh.pial.T1 ./tmp/rh_pial.gii
# mris_convert ./freesurfer/cv_01_full_brain/surf/lh.white ./tmp/lh_white.gii
# mris_convert ./freesurfer/cv_01_full_brain/surf/rh.white ./tmp/rh_white.gii
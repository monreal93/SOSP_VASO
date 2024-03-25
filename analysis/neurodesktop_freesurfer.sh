cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data


folder="05192023_sv_paper"
subject="sv_01"
scan="registered_Warped_sv"


cd ${folder}
cd ./analysis/
mkdir ./freesurfer

export SINGULARITYENV_SUBJECTS_DIR=./freesurfer

# Neurodesktop tools to use:
ml freesurfer

t1_file=./ants/${scan}.nii.gz

recon-all -autorecon1 -notal-check -skullstrip -wsthresh 5 -clean-bm -no-wsgcaatlas -subject ${subject} -i ${t1_file} #./dicom_mp2rage_iso0.7mm_iPAT3_20230519112423_17.nii

recon-all -autorecon2 -subject ${subject}

recon-all -autorecon3 -subject ${subject}

# If WM segmentation fails, add some control points and thren recon again...
# https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/ControlPoints_freeview/
# recon-all -autorecon2-cp -autorecon3 -subjid cp_before

# # saving WM and GM as masks..
mri_extract_label aseg.mgz 2 41 ../../../sv_01/wm_mask.nii
mri_extract_label aseg.mgz 3 42 ../../../sv_01/gm_mask.nii
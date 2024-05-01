cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="05252023_sv_paper"

scan="cv_41"
# scan_t1="dicom_mp2rage_iso0.7mm_iPAT3_20230519112423_17"
scan_t1="brain"

##### README.... 
# Before running this script:
# - do an intial registration in ITKsnap... save the values in a file called initial_matrix_$scan
# - Also create a ROI mask, that will be used for registration... save it as ./scan/scan_roi_msk.nii

cd ${folder}
cd ./analysis/
mkdir ./ants
chmod ugo+rwx ./ants/
cd ./ants

# Neurodesktop tools to use:
ml ants
ml afni

t1_file=../t1/${scan_t1}.nii
t1_vaso_file=../${scan}/T1_msk.nii # original
# t1_vaso_file=../${scan}/mean_v_msk.nii # original
mask=../${scan}/${scan}_roi_msk.nii     # Mask created from ITKsnap..
output1=registered_Warped_${scan}.nii
output2=registered_InverseWarped_${scan}.nii
initial_mtx=initial_matrix_${scan}.txt                     

# # Basic registration after registration in ITKsnap... save the values in a file called initial_matrix
# antsApplyTransforms --interpolation BSpline[5] -d 3 \
# -i "${t1_file}" \
# -r "${t1_vaso_file}" \
# -t "${initial_mtx}" \
# -o T1_VASO_registered_sample.nii

# # ANTs registration fMRI to Anatomical
# antsRegistration \
# --verbose 1 \
# --dimensionality 3 \
# --float 1 \
# --output [registered_,registered_Warped.nii.gz,registered_InverseWarped.nii.gz] \
# --interpolation Linear \
# --use-histogram-matching 0 \
# --winsorize-image-intensities [0.005,0.995] \
# --transform Similarity[0.05] \
# --initial-moving-transform initial_matrix.txt \
# --metric MI["${t1_file}","${t1_vaso_file}",0.7,32,Regular,0.1] \
# --convergence [1000x500,1e-6,10] \
# --shrink-factors 2x1 \
# --smoothing-sigmas 1x0vox \
# -x ./mask.nii

# --initial-moving-transform initial_matrix_t1_sv.txt \

# # ANTs registration Anatomical to spiral fMRI
# antsRegistration \
# --verbose 1 \
# --dimensionality 3 \
# --float 1 \
# --output [registered_,registered2_Warped.nii.gz,registered2_InverseWarped.nii.gz] \
# --interpolation Linear \
# --use-histogram-matching 0 \
# --winsorize-image-intensities [0.005,0.995] \
# --transform Similarity[0.05] \
# --metric MI["${t1_vaso_file}","${t1_file}",0.7,32,Regular,0.1] \
# --convergence [1000x500,1e-6,10] \
# --shrink-factors 2x1 \
# --smoothing-sigmas 1x0vox \
# -x ./mask_t1_cv.nii

### From Renzo's page: ORIGINAL
antsRegistration \
--verbose 1 \
--dimensionality 3 \
--float 1 \
--output ["${scan}"_,"${output1}","${output2}"] \
--interpolation Linear \
--use-histogram-matching 0 \
--winsorize-image-intensities [0.005,0.995] \
--initial-moving-transform "${initial_mtx}" \
--transform Rigid[0.05] \
--metric CC["${t1_vaso_file}","${t1_file}",0.7,32,Regular,0.1] \
--convergence [1000x500,1e-6,10] \
--shrink-factors 2x1 \
--smoothing-sigmas 1x0vox \
--transform Affine[0.1] \
--metric MI["${t1_vaso_file}","${t1_file}",0.7,32,Regular,0.1] \
--convergence [1000x500,1e-6,10] \
--shrink-factors 2x1 \
--smoothing-sigmas 1x0vox \
--transform SyN[0.1,2,0] \
--metric CC["${t1_vaso_file}","${t1_file}",1,2] \
--convergence [500x100,1e-6,10] \
--shrink-factors 2x1 \
--smoothing-sigmas 1x0vox \
-x "${mask}"


#####  Applying transformation to T1, GM and WM mask
cp ./${output1} ../${scan}/${scan}_t1.nii

transform=${scan}_0GenericAffine.mat   
# GM mask Aplying advanced transform obrained from antsRegistration
antsApplyTransforms --interpolation BSpline[5] -d 3 \
-i ../t1/gm_msk.nii \
-r "${t1_vaso_file}" \
-t "${transform}" \
-o "${scan}"_gm_msk_reg.nii

cp "${scan}"_gm_msk_reg.nii ../${scan}/${scan}_gm_msk.nii

# WM mask Aplying advanced transform obrained from antsRegistration
antsApplyTransforms --interpolation BSpline[5] -d 3 \
-i ../t1/wm_msk.nii \
-r "${t1_vaso_file}" \
-t "${transform}" \
-o "${scan}"_wm_msk_reg.nii

cp "${scan}"_wm_msk_reg.nii ../${scan}/${scan}_wm_msk.nii

# Remove warp from T1 registered to EPI sapce dataset...
# 3drefit -nowarp registered_Warped_cv.nii.gz

# ### From Renzo's page: ORIGINAL
# antsRegistration \
# --verbose 1 \
# --dimensionality 3 \
# --float 1 \
# --output [registered_,registered_Warped.nii.gz,registered_InverseWarped.nii.gz] \
# --interpolation Linear \
# --use-histogram-matching 0 \
# --winsorize-image-intensities [0.005,0.995] \
# --initial-moving-transform initial_matrix.txt \
# --transform Rigid[0.05] \
# --metric CC["${t1_vaso_file}","${t1_file}",0.7,32,Regular,0.1] \
# --convergence [1000x500,1e-6,10] \
# --shrink-factors 2x1 \
# --smoothing-sigmas 1x0vox \
# --transform Affine[0.1] \
# --metric MI["${t1_vaso_file}","${t1_file}",0.7,32,Regular,0.1] \
# --convergence [1000x500,1e-6,10] \
# --shrink-factors 2x1 \
# --smoothing-sigmas 1x0vox \
# --transform SyN[0.1,2,0] \
# --metric CC["${t1_vaso_file}","${t1_file}",1,2] \
# --convergence [500x100,1e-6,10] \
# --shrink-factors 2x1 \
# --smoothing-sigmas 1x0vox \
# -x mask.nii
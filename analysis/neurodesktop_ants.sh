cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="05192023_sv_paper"

scan="cv_01"
scan_t1="dicom_mp2rage_iso0.7mm_iPAT3_20230519112423_17"

##### Before running this script do a intial registration in ITKsnap... save the values in a file called initial_matrix_$scan

cd ${folder}
cd ./analysis/
mkdir ./ants
cd ./ants

# Neurodesktop tools to use:
ml ants
ml afni

t1_file=../../raw/nifti/${scan_t1}.nii
t1_vaso_file=../${scan}/T1_weighted_masked.nii
mask=../${scan}/mask_t1_${scan}.nii
output1=registered_Warped_${scan}.nii
output2=registered_InverseWarped_${scan}.nii
initial_mtx=initial_matrix_${scan}.txt
mask=mask_${scan}.txt                       # Mask created from ITKsnap..

# Basic registration after registration in ITKsnap... save the values in a file called initial_matrix
antsApplyTransforms --interpolation BSpline[5] -d 3 \
-i "${t1_file}" \
-r "${t1_vaso_file}" \
-t "${initial_mtx}" \
-o T1_VASO_registered_sample.nii

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
--output [registered_,"${output1}","${output2}"] \
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
# Neurodesktop tools to use:
ml afni

cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="05252023_sv_paper"
scan="sv_41"
# r_a_tr=8     # rest/activity TRs paper=sv(8), cv(7)
# tr=1.66      # volume TR paper=sv(1.66), cv(1.81)

# # Spiral reconstruction options
# traj="_nom" && cs="_cs" && b0="_b0" && k0="_k0" && rDORK="_rDORK"

cd ${folder}
mkdir ./analysis/${scan}
chmod ugo+rwx ./analysis/${scan}
cd ./analysis/${scan}

#######  ) Quantitive fMRI analysis
## 1) Manually align ananotmy to fMRI in ITK-SNAP, and create mask for registration in ANTs
## 2) Run neurodesktop_ants.sh to perform registration of T1 to fMRI space
## 3) Run neurodesktop_freesurfer.sh to perform segmentations...

# # Let's resample the gm and wm mask from freesurfer to fMRI space
# 3dresample -master ./T1_weighted_masked.nii -prefix gm_mask.nii -input gm_mask.nii -overwrite
# 3dresample -master ./T1_weighted_masked.nii -prefix wm_mask.nii -input wm_mask.nii -overwrite

# Let's make them binary mask, sometimes I use 0.5, sometimes 50, always check and redjust the "threshold..."
3dcalc -a ./"$scan"_gm_msk.nii -expr 'ispositive(a-50)' -prefix "$scan"_gm_msk1.nii -overwrite
3dcalc -a ./"$scan"_wm_msk.nii -expr 'ispositive(a-80)' -prefix "$scan"_wm_msk1.nii -overwrite

# Let's mask to a ROI in visual cortex
3dcalc -a "$scan"_gm_msk1.nii -b "$scan"_roi_msk.nii -expr 'a*b' -prefix "$scan"_gm_msk1.nii -overwrite
3dcalc -a "$scan"_wm_msk1.nii -b "$scan"_roi_msk.nii -expr 'a*b' -prefix "$scan"_wm_msk1.nii -overwrite

# Let's split the activations from GM and WM
# ToDo: Not sure if I want to use clustered_BOLD or BOLD_mask
3dcalc -a ./"$scan"_gm_msk1.nii -b ./BOLD_msk.nii -expr 'a*b' -prefix BOLD_gm.nii -overwrite
3dcalc -a ./"$scan"_wm_msk1.nii -b ./BOLD_msk.nii -expr 'a*b' -prefix BOLD_wm.nii -overwrite
3dcalc -a ./"$scan"_gm_msk1.nii -b ./VASO_msk.nii -expr 'a*b' -prefix VASO_gm.nii -overwrite
3dcalc -a ./"$scan"_wm_msk1.nii -b ./VASO_msk.nii -expr 'a*b' -prefix VASO_wm.nii -overwrite

# Let's ommit the first two slices (fold over artifacts)
slices=$(3dinfo -nk "$scan"_gm_msk.nii)
slices=$(($slices-1))
3dZcutup -keep 2 "${slices}" -prefix BOLD_gm.nii BOLD_gm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix BOLD_wm.nii BOLD_wm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix VASO_gm.nii VASO_gm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix VASO_wm.nii VASO_wm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix "$scan"_gm_msk1.nii "$scan"_gm_msk1.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix "$scan"_wm_msk1.nii "$scan"_wm_msk1.nii -overwrite

rm BOLD_roc*
rm VASO_roc*
step=0.3
##### WM VASO ROC
echo -e "Threshold \t Mean \t Voxels" >> VASO_roc_wm.txt
3dcalc -overwrite -a "$scan"_wm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_wm_msk1.nii)
echo -e "\n 99   $VASO_roc" >> VASO_roc_wm.txt

thr=0
for (( i=1; i<=20; i++))
do
    3dcalc -overwrite -a VASO_wm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet VASO_wm.nii)
    echo -e "\n $thr   $VASO_roc" >> VASO_roc_wm.txt
    thr=$(echo $thr+$step | bc -l)
done

##### GM VASO ROC
echo -e "Threshold \t Mean \t Voxels" >> VASO_roc_gm.txt
3dcalc -overwrite -a "$scan"_gm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_gm_msk1.nii)
echo -e "\n 99   $VASO_roc" >> VASO_roc_gm.txt

thr=0
for (( i=1; i<=20; i++))
do
    3dcalc -overwrite -a VASO_gm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet VASO_gm.nii)
    echo -e "\n $thr   $VASO_roc" >> VASO_roc_gm.txt
    thr=$(echo $thr+$step | bc -l)
done

##### WM BOLD ROC
echo -e "Threshold \t Mean \t Voxels" >> BOLD_roc_wm.txt
3dcalc -overwrite -a "$scan"_wm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_wm_msk1.nii)
echo -e "\n 99   $BOLD_roc" >> BOLD_roc_wm.txt

thr=0
for (( i=1; i<=20; i++))
do
    3dcalc -overwrite -a BOLD_wm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet BOLD_wm.nii)
    echo -e "\n $thr   $BOLD_roc" >> BOLD_roc_wm.txt
    thr=$(echo $thr+$step | bc -l)
done

##### GM BOLD ROC
echo -e "Threshold \t Mean \t Voxels" >> BOLD_roc_gm.txt
3dcalc -overwrite -a "$scan"_gm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_gm_msk1.nii)
echo -e "\n 99   $BOLD_roc" >> BOLD_roc_gm.txt

thr=0
for (( i=1; i<=20; i++))
do
    3dcalc -overwrite -a BOLD_gm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet BOLD_gm.nii)
    echo -e "\n $thr   $BOLD_roc" >> BOLD_roc_gm.txt
    thr=$(echo $thr+$step | bc -l)
done
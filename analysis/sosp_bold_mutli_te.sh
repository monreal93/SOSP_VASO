# Neurodesktop tools to use:
ml afni

cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="02132025_sb_9T_paper"
scan="sb_122_DS_SO_06mm_18fovz_12te_6te"
traj="girf"   # nom or girf
thr=1.6                    # z-score threshold

cd ./${folder}/analysis

######## Echo 1
cd ./${scan}_${traj}_ech1
### Percentage of active voxels
active_vox=$(3dBrickStat -count -non-zero -mask ./roi_mask.nii ./clustered_BOLD.nii)
mask_vox=$(3dBrickStat -count -non-zero -mask ./roi_mask.nii ./roi_mask.nii)
active_vox_per=$(echo $active_vox/$mask_vox*100 | bc -l)
echo -e "ROI Percentage active voxels threshold = $thr \n $active_vox_per" >> results.txt

### Median z-score
3dcalc -overwrite -a clustered_BOLD.nii -b roi_mask.nii -expr "a*b" -prefix clustered_BOLD_ROI.nii
med_z_score=$(3dBrickStat -median -non-zero -mask ./roi_mask.nii ./clustered_BOLD.nii)
echo -e "ROI Median z-score: \n $med_z_score" >> results.txt

####### ) Mean tSNR and effective tSNR
3dcalc -overwrite -a tSNR_msk.nii -b roi_mask.nii -expr "a*b" -prefix tSNR_ROI.nii
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./roi_mask.nii ./tSNR_msk.nii) 
echo -e "ROI BOLD brain mean tSNR: \n $mean_tSNR_b" >> results.txt


######## Echo 2
cd ../${scan}_${traj}_ech2
### Percentage of active voxels
active_vox=$(3dBrickStat -count -non-zero -mask ./roi_mask.nii ./clustered_BOLD.nii)
mask_vox=$(3dBrickStat -count -non-zero -mask ./roi_mask.nii ./roi_mask.nii)
active_vox_per=$(echo $active_vox/$mask_vox*100 | bc -l)
echo -e "ROI Percentage active voxels threshold = $thr \n $active_vox_per" >> results.txt

### Median z-score
3dcalc -overwrite -a clustered_BOLD.nii -b roi_mask.nii -expr "a*b" -prefix clustered_BOLD_ROI.nii
med_z_score=$(3dBrickStat -median -non-zero -mask ./roi_mask.nii ./clustered_BOLD.nii)
echo -e "ROI Median z-score: \n $med_z_score" >> results.txt

####### ) Mean tSNR and effective tSNR
3dcalc -overwrite -a tSNR_msk.nii -b roi_mask.nii -expr "a*b" -prefix tSNR_ROI.nii
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./roi_mask.nii ./tSNR_msk.nii) 
echo -e "ROI BOLD brain mean tSNR: \n $mean_tSNR_b" >> results.txt

###### overlapping btw different echos.... on ROI mask...
cd ..
mkdir ${scan:0:6}_te_comparison
cd ${scan:0:6}_te_comparison
3dcalc -overwrite -a ../${scan}_${traj}_ech1/clustered_BOLD.nii -b ../${scan}_${traj}_ech1/roi_mask.nii -expr "step(a-$thr)*b" -prefix te1_bin_output.nii
3dcalc -overwrite -a ../${scan}_${traj}_ech2/clustered_BOLD.nii -b ../${scan}_${traj}_ech1/roi_mask.nii -expr "step(a-$thr)*b" -prefix te2_bin_output.nii

3dcalc -overwrite -a te1_bin_output.nii -b te2_bin_output.nii -exp "a*b" -prefix te1_te2_bin_output.nii

3dcalc -overwrite -a te1_te2_bin_output.nii -b ../${scan}_${traj}_ech1/clustered_BOLD.nii -expr "a*b" -prefix BOLD_te1_common.nii
3dcalc -overwrite -a te1_te2_bin_output.nii -b ../${scan}_${traj}_ech2/clustered_BOLD.nii -expr "a*b" -prefix BOLD_te2_common.nii
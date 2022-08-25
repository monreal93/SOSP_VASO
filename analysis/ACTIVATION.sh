#!/bin/bash


echo "I estimate activity"

echo "based on difference" 
# AMM: Here I would need to replace the 9-19 and 29-39 with the actual volumes where I have activation...

# 3dTstat -mean -prefix VASO_r.nii -overwrite  VASO_trialAV_LN.nii'[9-19]'
# 3dTstat -mean -prefix VASO_a.nii -overwrite  VASO_trialAV_LN.nii'[29-39]'
3dTstat -mean -prefix VASO_r.nii -overwrite  VASO_trialAV_LN.nii'[5-11]'
3dTstat -mean -prefix VASO_a.nii -overwrite  VASO_trialAV_LN.nii'[12-17]'
3dcalc -a  VASO_r.nii -b VASO_a.nii -overwrite -expr '(a-b)/a' -prefix delta_VASO.nii

# Here  I need to change TR_times, and stim_times, I
echo "based on GLM" 
3dDeconvolve -overwrite -jobs 16 -polort a -input VASO_LN.nii\
             -num_stimts 1 \
             # -TR_times 3.6382 \
             # -stim_times 1 '1D: 0 60 120 180 240 300 360 420 480 540 600 660' 'UBLOCK(30,1)' -stim_label 1 Task \
             -TR_times 3.6382 \
             -stim_times 1 '1D: 0 60 120 180 240 300 360 420 480 540 600 660' 'UBLOCK(30,1)' -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_VASO.nii \
             -bucket STATS_VASO.nii

3dcalc -a HRF_VASO.nii'[1]'    -expr 'a'    -prefix 1_HRF_VASO.nii   -overwrite 
3dcalc -a STATS_VASO.nii'[0]'  -expr 'a'    -prefix 0_STATS_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[2]'  -expr 'a'    -prefix 2_STATS_VASO.nii -overwrite 


echo "The same for BOLD" 
3dTstat -mean -prefix BOLD_r.nii -overwrite  BOLD_trialAV_LN.nii'[5-11]'
3dTstat -mean -prefix BOLD_a.nii -overwrite  BOLD_trialAV_LN.nii'[17-23]'
3dcalc -a  BOLD_r.nii -b BOLD_a.nii -overwrite -expr '(b-a)/a' -prefix delta_BOLD.nii



3dDeconvolve -overwrite -jobs 16 -polort a -input BOLD_intemp.nii\
             -num_stimts 1 \
             -TR_times 1.5 \
             -stim_times 1 '1D: 0 60 120 180 240 300 360 420 480 540 600 660' 'UBLOCK(30,1)' -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_BOLD.nii \
             -bucket STATS_BOLD.nii

3dcalc -a HRF_BOLD.nii'[1]'    -expr '-1*a' -prefix 1_HRF_BOLD.nii -overwrite
3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'   -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_BOLD.nii -overwrite 

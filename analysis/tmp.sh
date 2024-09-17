


#### BOLD, Motion ############################
echo "BOLD based on GLM..."
3dDeconvolve -overwrite -jobs 16 -polort a \
            -force_TR $tr \
            -input ./${scan}_b_ups_mc_hpf.nii\
            -num_stimts 7 \
            -TR_times $tr \
            -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
            -tout \
            -x1D MODEL_wm \
            -iresp 1 HRF_BOLD.nii \
            -errts residual_BOLD.nii \
            -stim_file 2 motion_b.1D'[0]' -stim_base 2 -stim_label 2 roll \
            -stim_file 3 motion_b.1D'[1]' -stim_base 3 -stim_label 3 pitch \
            -stim_file 4 motion_b.1D'[2]' -stim_base 4 -stim_label 4 yaw \
            -stim_file 5 motion_b.1D'[3]' -stim_base 5 -stim_label 5 dS \
            -stim_file 6 motion_b.1D'[4]' -stim_base 6 -stim_label 6 dL \
            -stim_file 7 motion_b.1D'[5]' -stim_base 7 -stim_label 7 dP \
            -bucket STATS_BOLD.nii


######## BOLD
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a'    -prefix 1_HRF_BOLD.nii   -overwrite 
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_BOLD.nii   -overwrite 

3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'    -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_BOLD.nii -overwrite
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr 'a'    -prefix 2_STATS_BOLD.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii ./${scan}_b_ups_mc_hpf.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

# Masking
3dcalc -a 2_STATS_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii -overwrite


#### BOLD NO motion ########################
echo "BOLD based on GLM..."
3dDeconvolve -overwrite -jobs 16 -polort a \
            -force_TR $tr \
            -input ./${scan}_b_ups_mc_hpf.nii\
            -num_stimts 1 \
            -TR_times $tr \
            -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
            -tout \
            -x1D MODEL_wm \
            -iresp 1 HRF_BOLD.nii \
            -errts residual_BOLD.nii \
            -bucket STATS_BOLD.nii


######## BOLD
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a'    -prefix 1_HRF_BOLD.nii   -overwrite 
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_BOLD.nii   -overwrite 

3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'    -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_BOLD.nii -overwrite
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr 'a'    -prefix 2_STATS_BOLD.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii ./${scan}_b_ups_mc_hpf.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

# Masking
3dcalc -a 2_STATS_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii -overwrite
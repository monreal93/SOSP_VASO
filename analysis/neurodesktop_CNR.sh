# Neurodesktop tools to use:
ml afni


folder="02132025_sb_9T_paper"
scan="sb_122_DS_SO_06mm_18fovz_12te_6te"
suffix="girf_ech2"

block_tr=14   # Block length in TRs
blocks=11

cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data/${folder}/analysis/${scan}_${suffix}/

3dZeropad -master mean.nii -prefix b_bin_output.nii b_bin_output.nii -overwrite

# I need to get the timeseries mean per block... For now manually subsetting the timeseries
3dMean -prefix mean_block.nii -overwrite \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*0 | bc)'..'$(echo $block_tr*1-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*1 | bc)'..'$(echo $block_tr*2-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*2 | bc)'..'$(echo $block_tr*3-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*3 | bc)'..'$(echo $block_tr*4-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*4 | bc)'..'$(echo $block_tr*5-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*5 | bc)'..'$(echo $block_tr*6-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*6 | bc)'..'$(echo $block_tr*7-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*7 | bc)'..'$(echo $block_tr*8-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*8 | bc)'..'$(echo $block_tr*9-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*9 | bc)'..'$(echo $block_tr*10-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*10 | bc)'..'$(echo $block_tr*11-1 | bc)']' \
    ${scan}_mc_hpf.nii'['$(echo $block_tr*11 | bc)'..'$(echo $block_tr*12-1 | bc)']'

# Getting mean rest and activation
# 3dTstat -mean -prefix mean_rest.nii mean_block.nii'[0..'$(echo $block_tr/2-1 | bc)']' -overwrite
# 3dTstat -mean -prefix mean_activity.nii mean_block.nii'['$(echo $block_tr/2 | bc)'..$]' -overwrite

# Ommiting the first TRs of each block to avoid the initial transient
3dTstat -mean -prefix mean_rest.nii mean_block.nii'[2..'$(echo $block_tr/2-2 | bc)']' -overwrite
3dTstat -mean -prefix mean_activity.nii mean_block.nii'['$(echo $block_tr/2+1 | bc)'..'$(echo $block_tr-3 | bc)']' -overwrite

# Getting delta signal
3dcalc -a mean_rest.nii -b mean_activity.nii -expr 'b-a' -prefix del_sig.nii -overwrite
3dcalc -a del_sig.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix del_sig_msk.nii -overwrite

# Getting percent signal change
3dcalc -a del_sig.nii -b mean_rest.nii -expr '(a/b)*100' -prefix psc.nii -overwrite
3dcalc -a psc.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix psc_msk.nii -overwrite

# Getting standard deviation
3dTstat -nzstdev -prefix std.nii ${scan}_mc_hpf.nii -overwrite
3dcalc -a std.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix std_msk.nii -overwrite

# Getting delta signal over std
3dcalc -a del_sig.nii -b std.nii -expr 'a/b' -prefix del_sig_std.nii -overwrite
3dcalc -a del_sig_std.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix del_sig_std_msk.nii -overwrite

# Getting standard deviation of resiudal
3dTstat -nzstdev -prefix std_res_bold.nii residual_BOLD.nii -overwrite
3dcalc -a std_res_bold.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix std_res_bold_msk.nii -overwrite

# Getting delta signal over std
3dcalc -a del_sig.nii -b std_res_bold.nii -expr 'a/b' -prefix del_sig_std_res_bold.nii -overwrite
3dcalc -a del_sig_std_res_bold.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix del_sig_std_res_bold_msk.nii -overwrite


########### Temp... trying some stuff....###########
# Timeseries over std..
3dcalc -a ${scan}_mc_hpf.nii -b std_res_bold.nii -expr 'a/b' -prefix timeseries_std.nii -overwrite
3dcalc -a timeseries_std.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix timeseries_std_msk.nii -overwrite

# I need to get the timeseries mean per block... For now manually subsetting the timeseries
3dMean -prefix timeseries_std_mean_block.nii -overwrite \
    timeseries_std_msk.nii'['$(echo $block_tr*0 | bc)'..'$(echo $block_tr*1-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*1 | bc)'..'$(echo $block_tr*2-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*2 | bc)'..'$(echo $block_tr*3-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*3 | bc)'..'$(echo $block_tr*4-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*4 | bc)'..'$(echo $block_tr*5-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*5 | bc)'..'$(echo $block_tr*6-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*6 | bc)'..'$(echo $block_tr*7-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*7 | bc)'..'$(echo $block_tr*8-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*8 | bc)'..'$(echo $block_tr*9-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*9 | bc)'..'$(echo $block_tr*10-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*10 | bc)'..'$(echo $block_tr*11-1 | bc)']' \
    timeseries_std_msk.nii'['$(echo $block_tr*11 | bc)'..'$(echo $block_tr*12-1 | bc)']'

3dROIstats -mask bin_roi_output.nii -1DRformat -quiet ./timeseries_std_msk.nii > timeseries_std_msk.dat

3dROIstats -mask bin_roi_output.nii -1DRformat -quiet ${scan}_mc_hpf.nii > timeseries_msk.dat


#####################################################
3dcalc -a mean_block.nii -b std_res_bold.nii -expr 'a/b' -prefix mean_block_std_res_bold.nii -overwrite

# Ommiting the first TRs of each block to avoid the initial transient
3dTstat -mean -prefix mean_rest_block_std_res_bold.nii mean_block_std_res_bold.nii'[2..'$(echo $block_tr/2-2 | bc)']' -overwrite
3dTstat -mean -prefix mean_activity_block_std_res_bold.nii mean_block_std_res_bold.nii'['$(echo $block_tr/2+1 | bc)'..'$(echo $block_tr-3 | bc)']' -overwrite

# Getting delta signal
3dcalc -a mean_rest_block_std_res_bold.nii -b mean_activity_block_std_res_bold.nii -expr 'b-a' -prefix del_sig_block_std_res_bold.nii -overwrite
3dcalc -a del_sig_block_std_res_bold.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix del_sig_block_std_res_bold_msk.nii -overwrite

3dROIstats -mask bin_roi_output.nii -1DRformat -quiet ./mean_block_std_res_bold.nii > mean_block_std_res_bold.dat
######
# Neurodesktop tools to use:
ml afni


folder="01152025_sb_9T_paper"
scan="sb_011_DS_SO_06mm_18fovz_12te_6te"
suffix="girf_ech2"

block_tr=14   # Block length in TRs
blocks=12

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
    ${scan}_mc_hpf.nii'['$(echo $block_tr*10 | bc)'..'$(echo $block_tr*11-1 | bc)']'
    # ${scan}_mc_hpf.nii'['$(echo $block_tr*11 | bc)'..'$(echo $block_tr*12-1 | bc)']'

# Getting mean rest and activation
# 3dTstat -mean -prefix mean_rest.nii mean_block.nii'[0..'$(echo $block_tr/2-1 | bc)']' -overwrite
# 3dTstat -mean -prefix mean_activity.nii mean_block.nii'['$(echo $block_tr/2 | bc)'..$]' -overwrite

# Ommiting the first TRs of each block to avoid the initial transient
3dTstat -mean -prefix mean_rest.nii mean_block.nii'[3..'$(echo $block_tr/2-1 | bc)']' -overwrite
3dTstat -mean -prefix mean_activity.nii mean_block.nii'['$(echo $block_tr/2+3 | bc)'..$]' -overwrite

# Getting percent signal change
3dcalc -a mean_rest.nii -b mean_activity.nii -expr '((b-a)/a)*100' -prefix psc.nii -overwrite
3dcalc -a psc.nii -b roi_mask.nii -c b_bin_output.nii -expr 'a*b*c' -prefix psc_msk.nii -overwrite

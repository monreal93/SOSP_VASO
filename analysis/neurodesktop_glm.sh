# Neurodesktop tools to use:
ml afni

cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="06192024_sv_josh"
scan="sv_01"
r_a_tr=5     # rest/activity TRs paper=sv(8), cv(7)
tr=2.393    # volume TR paper=sv(1.66), cv(1.81)
motion_glm=0

# Spiral reconstruction options
traj="_nom" && cs="_cs" && b0="_b0" && k0="_k0" && rDORK="_rDORK"

cd ${folder}

mkdir ./analysis/${scan}
chmod ugo+rwx ./analysis/${scan}
cd ./analysis/${scan}

if [ "${scan:0:2}" = "sv" ]; then
    echo "Spiral VASO .."
    # Set path for reconstruction
    v_file=../../recon/${scan}_v${traj}${cs}${b0}${k0}${rDORK}.nii
    b_file=../../recon/${scan}_b${traj}${cs}${b0}${k0}${rDORK}.nii
    gre1=../../tmp/${scan}_1ech.nii
elif [ "${scan:0:2}" = "cv" ]; then
    echo "Cartesian VASO .."
    file=../../recon/${scan}_bv_epi.nii
    v_file=../../recon/${scan}_v_epi.nii
    b_file=../../recon/${scan}_b_epi.nii
fi

vol=$(3dinfo -nv ${v_file})

blocks=$(echo $vol/$r_a_tr/2 | bc -l)
blocks=$(echo ${blocks%.*})
block_dur=$(echo $tr*$blocks*2 | bc -l)
block_dur=$(echo ${block_dur%.*})
block_upsample=$(echo $blocks*2 | bc -l)
block_trs=$(echo $r_a_tr*2 | bc -l)
block_trs=$(echo ${block_trs%.*})

#### ) Activation maps
block_dur=$(echo $tr*$block_trs | bc -l)
tmp=0
start=1
# end=$(echo $blocks-1 | bc -l) # original
# stim_times='1D: 0 '  # original 
end=$(echo $blocks-2 | bc -l)  # new
stim_times=$(echo "1D: $block_dur ") # new
tmp=$(echo $tmp+$block_dur*2+$block_dur | bc -l)  # If temporal upsampling # New
stim_times=$(echo "$stim_times $tmp ") # new

for (( i=$start; i<=$end; i++))
do
	# tmp=$(echo $tmp+$block_dur*2+$block_dur | bc -l)  # If temporal upsampling # New
    tmp=$(echo $tmp+$block_dur*2 | bc -l)  # If temporal upsampling # Original
    # tmp=$(echo $tmp+$block_dur | bc -l)
	stim_times=$(echo "$stim_times $tmp ")
done
tmp=$(echo $block_dur | bc -l)
tmp=$(echo ${tmp%%.*})
ublock=$(echo "UBLOCK($tmp,1)")

if [ "$motion_glm" = 1 ]; then
    ### VASO GLM
    echo "VASO based on GLM..." 
    3dDeconvolve -overwrite -jobs 16 -polort a -input VASO_LN.nii\
                -num_stimts 7 \
                -TR_times $tr \
                -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
                -tout \
                -x1D MODEL_wm \
                -iresp 1 HRF_VASO.nii \
                -errts residual_VASO.nii \
                -stim_file 2 motion_v.1D'[0]' -stim_base 2 -stim_label 2 roll \
                -stim_file 3 motion_v.1D'[1]' -stim_base 3 -stim_label 3 pitch \
                -stim_file 4 motion_v.1D'[2]' -stim_base 4 -stim_label 4 yaw \
                -stim_file 5 motion_v.1D'[3]' -stim_base 5 -stim_label 5 dS \
                -stim_file 6 motion_v.1D'[4]' -stim_base 6 -stim_label 6 dL \
                -stim_file 7 motion_v.1D'[5]' -stim_base 7 -stim_label 7 dP \
                -bucket STATS_VASO.nii

    #### BOLD
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
else
    ### VASO GLM
    echo "VASO based on GLM..." 
    3dDeconvolve -overwrite -jobs 16 -polort a \
                -force_TR $tr \
                -input VASO_LN.nii\
                -num_stimts 1 \
                -TR_times $tr \
                -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
                -tout \
                -x1D MODEL_wm \
                -iresp 1 HRF_VASO.nii \
                -errts residual_VASO.nii \
                -bucket STATS_VASO.nii

    #### BOLD
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
fi

# VASO
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a'    -prefix 1_HRF_VASO.nii   -overwrite 
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_VASO.nii   -overwrite 

3dcalc -a STATS_VASO.nii'[0]'  -expr 'a'    -prefix 0_STATS_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_VASO.nii -overwrite
3dcalc -a STATS_VASO.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[2]'  -expr 'a'    -prefix 2_STATS_VASO.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii VASO_LN.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_VASO.nii      -expr 'b/a*100' -prefix 1_HRF_percent_VASO.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_VASO.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_VASO.nii

# BOLD
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
# # Original
# 3dcalc -a 2_STATS_VASO.nii -b mask.nii -expr 'a*b' -prefix VASO_msk.nii -overwrite
# 3dcalc -a 2_STATS_NEG_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii -overwrite
# New
3dcalc -a 2_STATS_NEG_VASO.nii -b mask.nii -expr 'a*b' -prefix VASO_msk.nii -overwrite
3dcalc -a 2_STATS_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii -overwrite

# Let's ommit the first two slices (fold over artifacts)
# ToDo: Find better way of doing this
slices=$(3dinfo -nk VASO_msk.nii)
slices=$(($slices-1))
# VASO
3dZcutup -keep 2 "${slices}" -prefix VASO_msk.nii VASO_msk.nii -overwrite
3dZeropad -I 2 -prefix VASO_msk.nii VASO_msk.nii -overwrite
# 3dZeropad -A 2 -prefix VASO_msk.nii VASO_mskt.nii -overwrite # For cv, sometimes I need to add the dim here
# BOLD
3dZcutup -keep 2 "${slices}" -prefix BOLD_msk.nii BOLD_msk.nii -overwrite
3dZeropad -I 2 -prefix BOLD_msk.nii BOLD_msk.nii -overwrite
# 3dZeropad -A 2 -prefix BOLD_msk.nii BOLD_msk.nii -overwrite # For cv, sometimes I need to add the dim here

##### ) Cluster activations
# Here, I need to play with the values after -1clip:
# -1clip threshold.. (~1.8), (1.5,1.2,270)
# rmm = cluster connection radius, larger value->remove small clusters
# vmul minimum cluster volume, smaller value->removes small clusters
3dclust -1noneg -overwrite -prefix clustered_VASO.nii -1clip 2 1.4 120 VASO_msk.nii
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 4 1.4 120 BOLD_msk.nii

### VASO mean residual
3dcalc -overwrite -a clustered_VASO.nii -expr 'step(a-1.8)' -prefix v_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask v_bin_output.nii -1DRformat -quiet residual_VASO.nii > residual_VASO.dat

### BOLD mean residual
3dcalc -overwrite -a clustered_BOLD.nii -expr 'step(a-1.8)' -prefix b_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask b_bin_output.nii -1DRformat -quiet ./residual_BOLD.nii > residual_BOLD.dat

# Let's ommit the first two slices (fold over artifacts)
# ToDo: Find better way of doing this
slices=$(3dinfo -nk VASO_msk.nii)
slices=$(($slices-1))
# VASO
3dZcutup -keep 2 "${slices}" -prefix VASO_msk.nii VASO_msk.nii -overwrite
3dZeropad -I 2 -prefix VASO_msk.nii VASO_msk.nii -overwrite
# 3dZeropad -A 2 -prefix VASO_msk.nii VASO_mskt.nii -overwrite # For cv, sometimes I need to add the dim here
# BOLD
3dZcutup -keep 2 "${slices}" -prefix BOLD_msk.nii BOLD_msk.nii -overwrite
3dZeropad -I 2 -prefix BOLD_msk.nii BOLD_msk.nii -overwrite
# 3dZeropad -A 2 -prefix BOLD_msk.nii BOLD_msk.nii -overwrite # For cv, sometimes I need to add the dim here


##### ) Get mean of activations and percentage change
### VASO
3dcalc -overwrite -a clustered_VASO.nii -expr 'step(a-1.8)' -prefix v_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask v_bin_output.nii -1DRformat -quiet VASO_LN.nii > timecourse_VASO.dat

### BOLD
3dcalc -overwrite -a clustered_BOLD.nii -expr 'step(a-1.8)' -prefix b_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask b_bin_output.nii -1DRformat -quiet ./${scan}_b_ups_mc_hpf.nii > timecourse_BOLD.dat

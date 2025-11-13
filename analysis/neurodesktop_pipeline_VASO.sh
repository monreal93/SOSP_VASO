# Neurodesktop tools to use:
ml afni
ml fsl
ml laynii

cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="05072025_sv_7T_m1"
scan="sv_012302_DS_SO_fs_pp"
r_a_tr=6     # rest/activity TRs paper=sv(8), cv(7)
tr=2.453   # volume TR paper=sv(1.66), cv(1.81)

# Spiral reconstruction options
traj="_nom" && cs="_cs" && b0="_tsb0" && co="" && k0="" && rDORK="_rDORK"
suffix=""

cd ${folder}

mkdir ./analysis/${scan}
chmod ugo+rwx ./analysis/${scan}
cd ./analysis/${scan}

if [ "${scan:0:2}" = "sv" ]; then
    echo "Spiral VASO .."
    # Set path for reconstruction
    v_file=../../recon/${scan}_v${traj}${cs}${b0}${co}${k0}${rDORK}${suffix}.nii
    b_file=../../recon/${scan}_b${traj}${cs}${b0}${co}${k0}${rDORK}${suffix}.nii
    gre1=../../tmp/${scan}${traj}.nii
elif [ "${scan:0:2}" = "cv" ]; then
    echo "Cartesian VASO .."
    file=../../recon/${scan}_bv_epi.nii
    v_file=../../recon/${scan}_v_epi.nii
    b_file=../../recon/${scan}_b_epi.nii
    # Split VASO and BOLD
    3dTcat -prefix ${b_file} ${file}'[0..$(2)]' -overwrite
    3dTcat -prefix ${v_file} ${file}'[1..$(2)]' -overwrite
fi

# ToDo: Get vol from recon nifti and write tr in nifti.. (when merging after recon)
# vol=160  # sv_01/02 =140, sv_03=220
# tr=1.66  # sv_01/02 = 1.87 , sv_03=0.93
# r_a_tr=8  # sv_01/02 = 7 , sv_03=11

vol=$(3dinfo -nv ${v_file})
# tr=$(3dinfo -tr ${v_file})
# tr=$(echo $tr*1000 | bc -l)

blocks=$(echo $vol/$r_a_tr/2 | bc -l)
blocks=$(echo ${blocks%.*})
block_dur=$(echo $tr*$blocks*2 | bc -l)
block_dur=$(echo ${block_dur%.*})
block_upsample=$(echo $blocks*2 | bc -l)
block_trs=$(echo $r_a_tr*2 | bc -l)
block_trs=$(echo ${block_trs%.*})

# ToDo: Find an easy way to realign Anatomy as fMRI...
# Center and deoblique datasets
3dinfo -obliquity ${b_file} >> obliquity.txt 
3drefit -xorigin cen -yorigin cen -zorigin cen ${v_file}
3drefit -xorigin cen -yorigin cen -zorigin cen ${b_file}
3drefit -xorigin cen -yorigin cen -zorigin cen ${gre1}
3dWarp -deoblique ${v_file}
3dWarp -deoblique ${b_file}
3dWarp -deoblique ${gre1}

# ) Making sure orientation is LPI
3drefit -orient LPI ${v_file}
3drefit -orient LPI ${b_file}
3drefit -orient LPI ${gre1}
3drefit -orient LPI ${v_file}
3drefit -orient LPI ${b_file}
3drefit -orient LPI ${gre1}

# # # Sometimes I need to Flip in y direction for CV
# 3dLRflip -X -prefix ${v_file} ${v_file} -overwrite
# 3dLRflip -X -prefix ${b_file} ${b_file} -overwrite

# 1) Creating Mask
#--- with FSL
bet $gre1 ./tmp -Z -f 0.4 -g 0 -n -m
mv tmp_mask.nii.gz mask.nii.gz
gzip -d mask.nii.gz
# #--- with AFNI
# if [ "${scan:0:2}" = "sv" ]; then
#     3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${gre1} -overwrite
# else
#     3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${v_file} -overwrite
# fi
3drefit -xdel $(3dinfo -adi ${b_file}) mask.nii
3drefit -ydel $(3dinfo -adj ${b_file}) mask.nii
3drefit -zdel $(3dinfo -adk ${b_file}) mask.nii

# Replacing non-steady state volumes..
# ToDo: How many volumes should I replace?
3dTcat -prefix ${b_file} ${b_file}'[3..5]' ${b_file}'[3..$]' -overwrite
3dTcat -prefix ${v_file} ${v_file}'[3..5]' ${v_file}'[3..$]' -overwrite

#### ) Temporal upsampling
3dUpsample -overwrite  -datum short -prefix ./${scan}_b_ups.nii -n 2 -input ./${b_file}
3dUpsample -overwrite  -datum short -prefix ./${scan}_v_ups.nii -n 2 -input ./${v_file}
# Updating TR times
3drefit -TR $tr ./${scan}_b_ups.nii
3drefit -TR $tr ./${scan}_v_ups.nii

# ) Motion correction
# ToDo: Check if mask can be used...
# ToDo: Use same function as Renzo's new scripts... 3dAllineate
echo "Motion correction... BOLD..."
# --- 3dvolreg
# 3dvolreg -base 4 -heptic -zpad 1 -overwrite -prefix ./${scan}_b_ups_mc.nii -1Dfile motion_b.1D ${scan}_b_ups.nii
# --- 3dAllineate
3dTstat -mean -overwrite -prefix ./n_ref.nii ${b_file}'[0..3]'  # create reference
3dAllineate -1Dmatrix_save  matrix.aff12.1D -1Dparam_save   param.aff12 -cost lpa \
    -prefix ./${scan}_b_ups_mc.nii -base ${output_dir} ./n_ref.nii -source ${scan}_b_ups.nii \
    -weight ./mask.nii -warp shift_rotate -final wsinc5

echo "Motion correction... VASO..."
# --- 3dvolreg
# 3dvolreg -base 4 -heptic -zpad 1 -overwrite -prefix ./${scan}_v_ups_mc.nii -1Dfile motion_v.1D ${scan}_v_ups.nii
# --- 3dAllineate
3dTstat -mean -overwrite -prefix ./n_ref.nii ${v_file}'[2..3]'  # create reference
3dAllineate -1Dmatrix_save  matrix.aff12.1D -1Dparam_save   param.aff12 -cost lpa \
    -prefix ./${scan}_v_ups_mc.nii -base ${output_dir} ./n_ref.nii -source ${scan}_v_ups.nii \
    -weight ./mask.nii -warp shift_rotate -final wsinc5

#### ) Temporal filtering
# This commands can be used for Tmporal filtering in AFNI 3dDetrend, 3dBandpass
echo "Temporal High-pass filter..."
bptf=$(echo $block_dur/$tr/2 | bc -l)
bptf=$(echo ${bptf%.*})
fslmaths ./${scan}_b_ups_mc.nii -Tmean tempmean
fslmaths ./${scan}_b_ups_mc.nii -bptf $bptf -1 -add tempmean ./${scan}_b_ups_mc_hpf.nii
fslchfiletype NIFTI ${scan}_b_ups_mc_hpf.nii
fslmaths ./${scan}_v_ups_mc.nii -Tmean tempmean
fslmaths ./${scan}_v_ups_mc.nii -bptf $bptf -1 -add tempmean ./${scan}_v_ups_mc_hpf.nii
fslchfiletype NIFTI ${scan}_v_ups_mc_hpf.nii

#### ) BOLD correction
LN_BOCO -Nulled ./${scan}_v_ups_mc_hpf.nii -BOLD ./${scan}_b_ups_mc_hpf.nii -trialBOCO $block_upsample

##### ) Calculating T1
echo "calculating T1 ..."
# Combining BOLD and VASO..
3dTcat -prefix combined.nii ./${scan}_b_ups_mc_hpf.nii ./${scan}_v_ups_mc_hpf.nii -overwrite
3dTstat -cvarinv -prefix T1_weighted.nii -overwrite combined.nii
rm combined.nii

#### ) Quality metrics
echo "calculating Mean and tSNR maps ..."
3dTstat -mean -prefix mean_b.nii ./${scan}_b_ups_mc_hpf.nii'[1..$]' -overwrite
3dTstat -mean -prefix mean_v.nii ./VASO_LN.nii'[1..$]' -overwrite
3dTstat  -overwrite -cvarinv  -prefix tSNR_b.nii ./${scan}_b_ups_mc_hpf.nii'[1..$]'
3dTstat  -overwrite -cvarinv  -prefix tSNR_v.nii VASO_LN.nii'[1..$]'

#### ) Activation maps
block_dur=$(echo $tr*$block_trs | bc -l)
tmp=0
start=1
# end=$(echo $blocks-1 | bc -l) # original
# stim_times='1D: 0 '  # original 
end=$(echo $blocks-2 | bc -l)
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

### Temp.. Manually writing the stim times...
# stim_times="1D: 15.9200 50.3600 80.3600 114.8000 149.2400 183.6800 218.1200 252.5600 287.0000 321.4400"
# stim_times="1D: 31.8400  100.7200  160.7200  229.6000  298.4800  367.3600"

# Finding rest and activity volumes
r1=0
r1=$(echo ${r1%.*})
r2=$(echo $blocks-1 | bc -l)
r2=$(echo ${r2%.*})
a1=$(echo $r2+1 | bc -l)
a1=$(echo ${a1%.*})
a2=$(echo $block_upsample-1 | bc -l)
a2=$(echo ${a2%.*})

tr_av_r=$(echo "[$r1-$r2]")
tr_av_a=$(echo "[$a1-$a2]")

#### VASO based on difference
echo "VASO based on difference..."
3dTstat -mean -prefix VASO_r.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_r"
3dTstat -mean -prefix VASO_a.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_a"
3dcalc -a  VASO_r.nii -b VASO_a.nii -overwrite -expr '(a-b)/a' -prefix delta_VASO.nii

### VASO GLM
echo "VASO based on GLM..." 
3dDeconvolve -overwrite -jobs 16 -polort a -input VASO_LN.nii\
             -num_stimts 1 \
             -TR_times $tr \
             -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_VASO.nii \
             -errts residual_VASO.nii \
             -bucket STATS_VASO.nii
             
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a'    -prefix 1_HRF_VASO.nii   -overwrite 
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_VASO.nii   -overwrite 

3dcalc -a STATS_VASO.nii'[0]'  -expr 'a'    -prefix 0_STATS_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_VASO.nii -overwrite
3dcalc -a STATS_VASO.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[2]'  -expr 'a'    -prefix 2_STATS_VASO.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii VASO_LN.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_VASO.nii      -expr 'b/a*100' -prefix 1_HRF_percent_VASO.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_VASO.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_VASO.nii

#### BOLD based on difference
echo "BOLD based on difference... " 
3dTstat -mean -prefix BOLD_r.nii -overwrite  BOLD_trialAV_LN.nii"$tr_av_r"
3dTstat -mean -prefix BOLD_a.nii -overwrite  BOLD_trialAV_LN.nii"$tr_av_a"
3dcalc -a  BOLD_r.nii -b BOLD_a.nii -overwrite -expr '(b-a)/a' -prefix delta_BOLD.nii

#### BOLD
echo "BOLD based on GLM..."
3dDeconvolve -overwrite -jobs 16 -polort a -input ./${scan}_b_ups_mc_hpf.nii\
             -num_stimts 1 \
             -TR_times $tr \
             -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_BOLD.nii \
             -errts residual_BOLD.nii \
             -bucket STATS_BOLD.nii

3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a'    -prefix 1_HRF_BOLD.nii   -overwrite 
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_BOLD.nii   -overwrite 

3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'    -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_BOLD.nii -overwrite
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr 'a'    -prefix 2_STATS_BOLD.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii ./${scan}_b_ups_mc_hpf.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

### VASO mean residual
3dcalc -overwrite -a clustered_VASO.nii -expr 'step(a-1.8)' -prefix v_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask v_bin_output.nii -1DRformat -quiet residual_VASO.nii > residual_VASO.dat

### BOLD mean residual
3dcalc -overwrite -a clustered_BOLD.nii -expr 'step(a-1.8)' -prefix b_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask b_bin_output.nii -1DRformat -quiet ./residual_BOLD.nii > residual_BOLD.dat

####### ) Masking relevant volumes
3dcalc -a T1_weighted.nii -b mask.nii -expr 'a*b' -prefix T1_msk.nii -overwrite
3dcalc -a mean_v.nii -b mask.nii -expr 'a*b' -prefix mean_v_msk.nii -overwrite
3dcalc -a mean_b.nii -b mask.nii -expr 'a*b' -prefix mean_b_msk.nii -overwrite
3dcalc -a 2_STATS_NEG_VASO.nii -b mask.nii -expr 'a*b' -prefix VASO_msk.nii -overwrite
3dcalc -a 2_STATS_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii -overwrite
3dcalc -a tSNR_v.nii -b mask.nii -expr 'a*b' -prefix tSNR_v_msk.nii -overwrite
3dcalc -a tSNR_b.nii -b mask.nii -expr 'a*b' -prefix tSNR_b_msk.nii -overwrite

3dcalc -a delta_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_delta_msk.nii -overwrite
3dcalc -a delta_VASO.nii -b mask.nii -expr 'a*b' -prefix VASO_delta_msk.nii -overwrite

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

##### ) Get mean of activations and percentage change
### VASO
3dcalc -overwrite -a clustered_VASO.nii -expr 'step(a-1.8)' -prefix v_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask v_bin_output.nii -1DRformat -quiet VASO_LN.nii > timecourse_VASO.dat

### BOLD
3dcalc -overwrite -a clustered_BOLD.nii -expr 'step(a-1.8)' -prefix b_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask b_bin_output.nii -1DRformat -quiet ./${scan}_b_ups_mc_hpf.nii > timecourse_BOLD.dat

####### ) Mean tSNR and effective tSNR
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./tSNR_b_msk.nii) 
echo -e "BOLD brain mean tSNR: \n $mean_tSNR_b" >> results.txt
mean_tSNR_v=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./tSNR_v_msk.nii) 
echo -e "VASO brain mean tSNR: \n $mean_tSNR_v" >> results.txt

3dcalc -a ./tSNR_b_msk.nii -expr "a/sqrt($tr)" -prefix ./eff_tSNR_b.nii -overwrite
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./eff_tSNR_b.nii) 
echo -e "BOLD brain mean effective tSNR: \n $mean_tSNR_b" >> results.txt
3dcalc -a ./tSNR_v_msk.nii -expr "a/sqrt($tr)" -prefix ./eff_tSNR_v.nii -overwrite
mean_tSNR_v=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./eff_tSNR_v.nii) 
echo -e "VASO brain mean effective tSNR: \n $mean_tSNR_v" >> results.txt



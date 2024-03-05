cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="02022024_sv_7Ti"
scan="sv_01_35G_155SR"

# ToDo: Get vol from recon nifti and write tr in nifti.. (when merging after recon)
vol=140  # sv_01/02 =140, sv_03=220
tr=1.87  # sv_01/02 = 1.87 , sv_03=0.93
r_a_tr=7  # sv_01/02 = 7 , sv_03=11

# Reconstruction options
traj="_sk" && cs="_cs" && b0="_b0" && k0="_k0" && rDORK="_rDORK"

cd ${folder}
mkdir ./analysis/${scan}
cd ./analysis/${scan}

# Neurodesktop tools to use:
ml afni
ml fsl
ml laynii

blocks=$(echo $vol/$r_a_tr/2 | bc -l)
blocks=$(echo ${blocks%.*})
block_dur=$(echo $tr*$blocks*2 | bc -l)
block_dur=$(echo ${block_dur%.*})
block_upsample=$(echo $blocks*2 | bc -l)
block_trs=$(echo $r_a_tr*2 | bc -l)
block_trs=$(echo ${block_trs%.*})

# Set path for reconstruction
v_file=../../recon/${scan}_v${traj}${cs}${b0}${k0}${rDORK}.nii
b_file=../../recon/${scan}_b${traj}${cs}${b0}${k0}${rDORK}.nii

# 1) Skull strip
gre1=../../tmp/${scan}_1ech.nii
bet $gre1 ./tmp -Z -f 0.4 -g 0 -n -m
mv tmp_mask.nii.gz mask.nii.gz
gzip -d mask.nii.gz

# Replacing non-steady state volumes..
3dTcat -prefix ${b_file} ${b_file}'[4..7]' ${b_file}'[4..$]' -overwrite
3dTcat -prefix ${v_file} ${v_file}'[4..7]' ${v_file}'[4..$]' -overwrite

# ) Motion correction
# ToDo: Check if mask can be used...
echo "Motion correction... BOLD..."
3dvolreg -base 2 -heptic -zpad 1 -overwrite -prefix ./${scan}_b_mc.nii -1Dfile motion_b.1D ${b_file}
echo "Motion correction... VASO..."
3dvolreg -base 2 -heptic -zpad 1 -overwrite -prefix ./${scan}_v_mc.nii -1Dfile motion_v.1D ${v_file}

#### ) Temporal upsampling
3dUpsample -overwrite  -datum short -prefix ./${scan}_b_mc_ups.nii -n 2 -input ./${scan}_b_mc.nii
3dUpsample -overwrite  -datum short -prefix ./${scan}_v_mc_ups.nii -n 2 -input ./${scan}_v_mc.nii
# Updating TR times
3drefit -TR $tr ./${scan}_b_mc_ups.nii
3drefit -TR $tr ./${scan}_v_mc_ups.nii

#### ) Temporal filtering
echo "Temporal High-pass filter..."
bptf=$(echo $block_dur/$tr/2 | bc -l)
bptf=$(echo ${bptf%.*})
fslmaths ./${scan}_b_mc_ups.nii -Tmean tempmean
fslmaths ./${scan}_b_mc_ups.nii -bptf $bptf -1 -add tempmean ./${scan}_b_mc_ups_hpf.nii
fslchfiletype NIFTI ${scan}_b_mc_ups_hpf.nii
fslmaths ./${scan}_v_mc_ups.nii -Tmean tempmean
fslmaths ./${scan}_v_mc_ups.nii -bptf $bptf -1 -add tempmean ./${scan}_v_mc_ups_hpf.nii
fslchfiletype NIFTI ${scan}_v_mc_ups_hpf.nii

#### ) BOLD correction
LN_BOCO -Nulled ./${scan}_v_mc_ups_hpf.nii -BOLD ./${scan}_b_mc_ups_hpf.nii -trialBOCO $block_trs

#### ) Calculating T1
echo "calculating T1 ..."
NumVol=`3dinfo -nv ./${scan}_b_mc_ups_hpf.nii`
# 3dcalc -a ./${scan}_b_mc_ups_hpf.nii'[3..'`expr $NumVol - 2`']' -b  ./${scan}_v_mc_ups_hpf.nii'[3..'`expr $NumVol - 2`']' -expr 'a+b' -prefix combined.nii -overwrite
# ToDo: Optimize combination of vaso and bold
3dTcat -overwrite ./${scan}_b_mc_ups_hpf.nii'[0]' -prefix combined.nii
3dTcat -overwrite ./combined.nii ./${scan}_v_mc_ups_hpf.nii'[0]' -prefix combined.nii
# for (( i = 1; i <= ${NumVol}; i++ ))
# do
#     3dTcat -overwrite ./combined.nii ./${scan}_b_mc_ups_hpf.nii"[$i]" -prefix combined.nii
#     3dTcat -overwrite ./combined.nii ./${scan}_v_mc_ups_hpf.nii"[$i]" -prefix combined.nii
done
3dTstat -cvarinv -prefix T1_weighted.nii -overwrite combined.nii 
rm combined.nii

#### ) Quality metrics
echo "calculating Mean and tSNR maps ..."
3dTstat -mean -prefix mean_b.nii ./${scan}_b_mc_ups_hpf.nii'[1..$]' -overwrite
3dTstat -mean -prefix mean_v.nii ./VASO_LN.nii'[1..$]' -overwrite
3dTstat  -overwrite -cvarinv  -prefix tSNR_b.nii ./${scan}_b_mc_ups_hpf.nii'[1..$]'
3dTstat  -overwrite -cvarinv  -prefix tSNR_v.nii VASO_LN.nii'[1..$]'

#### ) Activation maps
block_dur=$(echo $tr*$block_trs | bc -l)
tmp=0
start=1
end=$(echo $blocks-1 | bc -l)
stim_times='1D: 0 '

for (( i=$start; i<=$end; i++))
do
	tmp=$(echo $tmp+$block_dur*2 | bc -l)  # If temporal upsampling
    # tmp=$(echo $tmp+$block_dur | bc -l)
	stim_times=$(echo "$stim_times $tmp ")
done
tmp=$(echo $block_dur | bc -l)
tmp=$(echo ${tmp%%.*})
ublock=$(echo "UBLOCK($tmp,1)")

# Finding rest and activity volumes
tmp=$(echo $block_trs/4 | bc -l)
tmp=$(echo ${tmp%.*})
r1=$(echo $tmp-1 | bc -l)
r1=$(echo ${r1%.*})
r2=$(echo $r1+$tmp | bc -l)
r2=$(echo ${r2%.*})
a1=$(echo $r2+$tmp | bc -l)
a1=$(echo ${a1%.*})
a2=$(echo $a1+$tmp | bc -l)
a2=$(echo ${a2%.*})
tr_av_r=$(echo "[$r1-$r2]")
tr_av_a=$(echo "[$a1-$a2]")

#### VASO based on difference
echo "VASO based on difference..."
3dTstat -mean -prefix VASO_r.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_r"
3dTstat -mean -prefix VASO_a.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_a"
3dcalc -a  VASO_r.nii -b VASO_a.nii -overwrite -expr '(a-b)/a' -prefix delta_VASO.nii

#### VASO GLM
echo "VASO based on GLM..." 
3dDeconvolve -overwrite -jobs 16 -polort a -input VASO_LN.nii\
             -num_stimts 1 \
             -TR_times $tr \
             -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_VASO.nii \
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
3dDeconvolve -overwrite -jobs 16 -polort a -input ./${scan}_b_mc_ups_hpf.nii\
             -num_stimts 1 \
             -TR_times $tr \
             -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_BOLD.nii \
             -bucket STATS_BOLD.nii

3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a'    -prefix 1_HRF_BOLD.nii   -overwrite 
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_BOLD.nii   -overwrite 

3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'    -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_BOLD.nii -overwrite
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr 'a'    -prefix 2_STATS_BOLD.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii /${scan}_b_mc_ups_hpf.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

##### ) Cluster activations
# Here, I need to play with the values after -1clip:
# -1clip threshold.. (~1.8), (1.5,1.2,270)
# rmm = cluster connection radius, larger value->remove small clusters
# vmul minimum cluster volume, smaller value->removes small clusters
3dclust -1noneg -overwrite -prefix clustered_VASO.nii -1clip 2 1.4 120 2_STATS_VASO.nii
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 3 1.4 120 2_STATS_NEG_BOLD.nii

3dclust -1noneg -overwrite -prefix clustered_VASO.nii -1clip 0.02 1 120 delta_VASO.nii
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 0.02 1 120 delta_BOLD.nii

####### 10) Masking relevant volumes
fslmaths T1_weighted.nii -mul mask.nii T1_msk.nii
fslmaths clustered_BOLD.nii -mul mask.nii BOLD_clust_msk.nii
fslmaths clustered_VASO.nii -mul mask.nii VASO_clust_msk.nii

fslmaths delta_BOLD.nii -mul mask.nii BOLD_delta_msk.nii
fslmaths delta_VASO.nii -mul mask.nii VASO_delta_msk.nii
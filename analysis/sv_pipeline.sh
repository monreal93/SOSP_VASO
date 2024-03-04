
# Make sure you change the directory to the one where you are running the analysis....
# sosp_vaso/data/"date_of_scan"/analysis/...

# If I want to combine 2 diff runs, I can do it with fslmaths, after moco I can average them, example:
# fslmaths ./Documents/PhD/PhD_2022/sosp_vaso/data/_10282022_sv/analysis/sv_02/Nulled_Basis_b.nii -add ./Documents/PhD/PhD_2022/sosp_vaso/data/_10282022_sv/analysis/sv_03/Nulled_Basis_b.nii -div 2 ./Documents/PhD/PhD_2022/sosp_vaso/data/_10282022_sv/analysis/sv_02_03/Nulled_Basis_b.nii 

#### 1) Let's first copy the files I need into this folder
cp /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/analysis/youtube/gnuplot_moco.txt ./gnuplot_moco.txt
cp /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/analysis/youtube/gnuplot_layers.txt ./gnuplot_layers.txt

# I am not sure if I also want to copy mocobatch.m
# cp /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/analysis/youtube/

#### 2) MOTION CORRECTION (modify cp path in mocobatch.m):
file=../*sv_02_bv*
vol=140  # s=160 / c=140
tr=1.87  # s=1.66 / c=1.81
r_a_tr=7  # TRs for rest and activity s=8 / c=7

blocks=$(echo $vol/$r_a_tr/2 | bc -l)
blocks=$(echo ${blocks%.*})
block_dur=$(echo $tr*$blocks*2 | bc -l)
block_dur=$(echo ${block_dur%.*})
block_upsample=$(echo $blocks*2 | bc -l)
block_trs=$(echo $r_a_tr*2 | bc -l)
block_trs=$(echo ${block_trs%.*})

path=$(ls -1 $file | head -n 1)
dir=$(pwd)
beast_path="$/mnt/dabeast/5T4/Alejandro/sosp_vaso/data/${dir:53:11}"

# Different ways of masking/extracting brain
# 3dSkullStrip -push_to_edge -no_avoid_eyes -mask_vol -prefix mask.nii -overwrite -input /mnt/dabeast/5T3/Alejandro/sosp_vaso/data/05192023_sv/tmp/sv_01_gre_1echo.nii
# 3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${path}

if [ ${file:4:2} == "sv" ] then
	gre1=${beast_path}/tmp/${file:4:5}_gre_1echo
	bet ${gre1:1:100} ./tmp -Z -f 0.4 -g 0 -n -m
	mv tmp_mask.nii.gz mask.nii.gz
	gzip -d mask.nii.gz
else
	3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${path}	
fi

cp $file.nii ./Basis_a.nii
3dTcat -prefix Basis_a.nii Basis_a.nii'[4..7]' Basis_a.nii'[4..$]' -overwrite
cp ./Basis_a.nii ./Basis_b.nii

3dinfo -nt Basis_a.nii >> NT.txt
3dinfo -nt Basis_b.nii >> NT.txt

###### 3) Motion Correction
## for non-vaso data, use script fn_mocobatch_flex1
matlab -nodesktop -nosplash -r "fn_mocobatch_flex "${dir}

##
gnuplot "gnuplot_moco.txt"

rm ./Basis_*

################### Start here if I combine 2 runs after moco with sv_start_moco_flex.sh
################### I can combine the two runs with fslmaths, example above

#### 4) BOLD CORRECTION (modify -trialBOCO, and 3drefit -TR):
3dcalc -a Nulled_Basis_b.nii'[1..$(2)]' -expr 'a' -prefix Nulled.nii -overwrite
3dcalc -a Not_Nulled_Basis_a.nii'[0..$(2)]' -expr 'a' -prefix BOLD.nii -overwrite
3dUpsample -overwrite  -datum short -prefix Nulled_intemp.nii -n 2 -input Nulled.nii
3dUpsample -overwrite  -datum short -prefix BOLD_intemp.nii   -n 2 -input   BOLD.nii
NumVol=`3dinfo -nv BOLD_intemp.nii`
3dTcat -overwrite -prefix Nulled_intemp.nii Nulled_intemp.nii'[0]' Nulled_intemp.nii'[0..'`expr $NumVol - 2`']' 

### 4.1) Temporal filtering
echo "Temporal High-pass filter"
bptf=$(echo $block_dur/$tr/2 | bc -l)
bptf=$(echo ${bptf%.*})
fslmaths BOLD_intemp.nii -Tmean tempmean
fslmaths BOLD_intemp.nii -bptf $bptf -1 -add tempmean BOLD_intemp
fslchfiletype NIFTI BOLD_intemp.nii.gz
fslmaths Nulled_intemp.nii -Tmean tempmean
fslmaths Nulled_intemp.nii -bptf $bptf -1 -add tempmean Nulled_intemp

rm BOLD_intemp.nii Nulled_intemp.nii

#### 4.2) Bold correction
/home/laynii/LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO $block_trs

## 
echo "I am correcting for the proper TR in the header"
3drefit -TR $tr BOLD_intemp.nii
3drefit -TR $tr VASO_LN.nii

echo "calculating T1 in EPI space"
NumVol=`3dinfo -nv Nulled_Basis_b.nii`
3dcalc -a Nulled_Basis_b.nii'[3..'`expr $NumVol - 2`']' -b  Not_Nulled_Basis_a.nii'[3..'`expr $NumVol - 2`']' -expr 'a+b' -prefix combined.nii -overwrite
3dTstat -cvarinv -prefix T1_weighted.nii -overwrite combined.nii 
rm combined.nii

##
echo "calculating Mean and tSNR maps"
3dTstat -mean -prefix mean_nulled.nii Nulled.nii -overwrite
3dTstat -mean -prefix mean_notnulled.nii BOLD.nii -overwrite
  3dTstat  -overwrite -mean  -prefix BOLD.Mean.nii \
     BOLD_intemp.nii'[1..$]'
  3dTstat  -overwrite -cvarinv  -prefix BOLD.tSNR.nii \
     BOLD_intemp.nii'[1..$]'
  3dTstat  -overwrite -mean  -prefix VASO.Mean.nii \
     VASO_LN.nii'[1..$]'
  3dTstat  -overwrite -cvarinv  -prefix VASO.tSNR.nii \
     VASO_LN.nii'[1..$]'

echo "curtosis and skew"
/home/laynii/LN_SKEW -input BOLD.nii
/home/laynii/LN_SKEW -input VASO_LN.nii

#### 5) ACTIVTION (modify -TR_times and -stim_times in VASO and BOLD...):
# TR_times = VolumeTR(can get it from params.gen.volTR)
# stim_times = time of TR times duration of TRs in BLOCK times 2(Vaso-Bold), then just double this till number of Blocks
# stim_times = 0 TR_times*(#ON_TRs+#OFF_TRs)*2
# UBLOCK = (TR_times*#ON_TRs*2,1)
block_dur=$(echo $tr*($block_trs) | bc -l)
tmp=0
start=1
end=$(echo $blocks-1 | bc -l)
stim_times='1D: 0 '

for (( i=$start; i<=$end; i++))
do
	tmp=$(echo $tmp+($block_dur*2) | bc -l)
	stim_times=$(echo "$stim_times $tmp ")
done
tmp=$(echo $block_dur | bc -l)
tmp=$(echo ${tmp%%.*})
ublock=$(echo "UBLOCK($tmp,1)")

echo "VASO based on difference"
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

# AMM: Here I would need to replace the 9-19 and 29-39 with the actual volumes 
3dTstat -mean -prefix VASO_r.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_r"
3dTstat -mean -prefix VASO_a.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_a"
3dcalc -a  VASO_r.nii -b VASO_a.nii -overwrite -expr '(a-b)/a' -prefix delta_VASO.nii

# For now, I still need to copy and paste the stim_times, but should make it automatic
# printf "copy and paste this \n $stim_times \n $ublock \n"

###
echo "VASO based on GLM" 
3dDeconvolve -overwrite -jobs 16 -polort a -input VASO_LN.nii\
             -num_stimts 1 \
             -TR_times $tr \
             -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_VASO.nii \
             -bucket STATS_VASO.nii

             #-stim_times 1 '1D: 0 53.12 106.24 159.36 212.48 265.6 318.72 371.84' 'UBLOCK(26.56,1)' -stim_label 1 Task \
             
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a'    -prefix 1_HRF_VASO.nii   -overwrite 
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_VASO.nii   -overwrite 

3dcalc -a STATS_VASO.nii'[0]'  -expr 'a'    -prefix 0_STATS_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_VASO.nii -overwrite
3dcalc -a STATS_VASO.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[2]'  -expr 'a'    -prefix 2_STATS_VASO.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii VASO_LN.nii

3dcalc -a mean.nii -b 1_HRF_VASO.nii      -expr 'b/a*100' -prefix 1_HRF_percent_VASO.nii
3dcalc -a mean.nii -b  1_HRF_NEG_VASO.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_VASO.nii

echo "BOLD based on difference" 
# AMM: Here I would need to replace the 9-19 and 29-39 with the actual volumes 
3dTstat -mean -prefix BOLD_r.nii -overwrite  BOLD_trialAV_LN.nii"$tr_av_r"
3dTstat -mean -prefix BOLD_a.nii -overwrite  BOLD_trialAV_LN.nii"$tr_av_a"
3dcalc -a  BOLD_r.nii -b BOLD_a.nii -overwrite -expr '(b-a)/a' -prefix delta_BOLD.nii

##
echo "BOLD based on GLM"
3dDeconvolve -overwrite -jobs 16 -polort a -input BOLD_intemp.nii\
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

3dTstat -mean -overwrite -prefix mean.nii BOLD_intemp.nii

3dcalc -a mean.nii -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

##### 6) Cluster activations
# Here, I need to play with the values after -1clip:
# -1clip threshold.. (~1.8), (1.5,1.2,270)
# rmm = cluster connection radius, larger value->remove small clusters
# vmul minimum cluster volume, smaller value->removes small clusters
3dclust -1noneg -overwrite -prefix clustered_VASO.nii -1clip 2 1.4 120 2_STATS_VASO.nii
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 3 1.4 120 2_STATS_NEG_BOLD.nii

####### 10) Masking relevant volumes
fslmaths T1_weighted.nii -mul mask.nii T1_weighted_masked.nii
fslmaths clustered_BOLD.nii -mul mask.nii clustered_BOLD_masked.nii
fslmaths clustered_VASO.nii -mul mask.nii clustered_VASO_masked.nii

### Create the manual Visual and Motor mask....
##### 11) Masking again for a ROI, if needed....
fslmaths mask.nii -mul roi_visual_motor.nii.gz roi_visual_motor_new.nii
fslmaths clustered_BOLD_masked.nii -mul roi_visual_motor_new.nii clustered_BOLD_masked.nii
fslmaths clustered_VASO_masked.nii -mul roi_visual_motor_new.nii clustered_VASO_masked.nii


##### 7) Get mean of activations and percentage change
### VASO
3dcalc -overwrite -a clustered_VASO_masked.nii -expr 'step(a-1.8)' -prefix binarised_output.nii
3dROIstats -mask binarised_output.nii -1DRformat -quiet VASO_LN.nii > timecourse_VASO.dat
# Percentage change (Modify the clip level in step(b-cliplevel)
# 3dTstat -prefix -overwrite mean_tmp.nii VASO_LN.nii
# 3dClipLevel mean_tmp.nii
# 3dcalc -a VASO_LN.nii -b mean_tmp.nii -expr "((((a-b)/1000) * 100) * step(b-0.345117))" -prefix -overwrite scaled.nii
# 3dROIstats -mask binarised_output.nii -1DRformat -quiet scaled.nii > vaso_per.dat

### BOLD
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'step(a-1.8)' -prefix binarised_output.nii
3dROIstats -mask binarised_output.nii -1DRformat -quiet BOLD_intemp.nii.gz > timecourse_BOLD.dat
# Percentage change (Modify the clip level in step(b-cliplevel)
# 3dTstat -prefix -overwrite mean_tmp.nii BOLD_intemp.nii.gz
# 3dClipLevel mean_tmp.nii
# 3dcalc -a BOLD_intemp.nii.gz -b mean_tmp.nii -expr "((((a-b)/1000) * 100) * step(b-180.7))" -prefix -overwrite scaled.nii
# 3dROIstats -mask binarised_output.nii -1DRformat -quiet scaled.nii > bold_per.dat

#### 8) Possible smoothing
LN_GRADSMOOTH -input 2_STATS_VASO.nii -gradfile T1_weighted.nii -FWHM 1.5 -within -selectivity 0.08 

##### 9) Layering
# Upsample
delta_x=$(3dinfo -di T1_weighted_masked.nii.gz)
delta_y=$(3dinfo -dj T1_weighted_masked.nii.gz)
delta_z=$(3dinfo -dk T1_weighted_masked.nii.gz)
sdelta_x=$(echo "((sqrt($delta_x * $delta_x) / 5))"|bc -l)
sdelta_y=$(echo "((sqrt($delta_y * $delta_y) / 5))"|bc -l)
sdelta_z=$(echo "((sqrt($delta_z * $delta_z) / 1))"|bc -l) 

## Scaling with AFNI... Not the best for oblique datasets
# here I only upscale in 2 dimensions. 
#3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode NN -overwrite -prefix scaled_$1 -input $1 
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode Cu -overwrite -prefix scaled_T1.nii -input T1_weighted_masked.nii.gz 
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode Cu -overwrite -prefix scaled_VASO.nii -input clustered_VASO_masked.nii.gz 
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode Cu -overwrite -prefix scaled_BOLD.nii -input clustered_BOLD_masked.nii.gz 


### NOW CREATE Manually the mask in FSLEYES, call it rim.nii

#estimating layers based on rim
LN_GROW_LAYERS -rim rim.nii -N 11 -vinc 40

#extractiong profiles for VASO
#get mean value, STDEV, and number of voxels
3dROIstats -mask rim_layers.nii -1DRformat -quiet -nzmean scaled_VASO.nii > layer_t.dat
3dROIstats -mask rim_layers.nii -1DRformat -quiet -sigma scaled_VASO.nii >> layer_t.dat
3dROIstats -mask rim_layers.nii -1DRformat -quiet -nzvoxels scaled_VASO.nii >> layer_t.dat
#format file to be in columns, so gnuplot can read it.
WRD=$(head -n 1 layer_t.dat|wc -w); for((i=2;i<=$WRD;i=i+2)); do awk '{print $'$i'}' layer_t.dat| tr '\n' ' ';echo; done > layer.dat

1dplot -sepscl layer.dat 
mv layer.dat layer_VASO.dat

#extractiong profiles for BOLD
#get mean value, STDEV, and number of voxels
3dROIstats -mask rim_layers.nii -1DRformat -quiet -nzmean scaled_BOLD.nii > layer_t.dat
3dROIstats -mask rim_layers.nii -1DRformat -quiet -sigma scaled_BOLD.nii >> layer_t.dat
3dROIstats -mask rim_layers.nii -1DRformat -quiet -nzvoxels scaled_BOLD.nii >> layer_t.dat
#format file to be in columns, so gnuplot can read it.
WRD=$(head -n 1 layer_t.dat|wc -w); for((i=2;i<=$WRD;i=i+2)); do awk '{print $'$i'}' layer_t.dat| tr '\n' ' ';echo; done > layer.dat

1dplot -sepscl layer.dat 
mv layer.dat layer_BOLD.dat

gnuplot "gnuplot_layers.txt"



######### FOR figures: Get activation ranges...


# VASO
3dcalc -overwrite -a clustered_VASO_masked.nii -expr 'within(a,2,3)' -prefix VASO_2_3.nii
3dcalc -overwrite -a clustered_VASO_masked.nii -expr 'within(a,3.001,4)' -prefix VASO_3_4.nii
3dcalc -overwrite -a clustered_VASO_masked.nii -expr 'within(a,4.001,5)' -prefix VASO_4_5.nii
3dcalc -overwrite -a clustered_VASO_masked.nii -expr 'within(a,5.001,6)' -prefix VASO_5_6.nii
3dcalc -overwrite -a clustered_VASO_masked.nii -expr 'within(a,6.001,50)' -prefix VASO_6_plus.nii

# BOLD
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,3,5)' -prefix BOLD_3_5.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,5.001,7)' -prefix BOLD_5_7.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,7.001,9)' -prefix BOLD_7_9.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,9.001,11)' -prefix BOLD_9_11.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,11.001,50)' -prefix BOLD_11_plus.nii


# Calculate volume mean tSNR
3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./BOLD_tSNR.nii
3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./VASO_LN_tSNR.nii






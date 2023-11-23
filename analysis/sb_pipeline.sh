
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
file=../*abc_01*
vol=180  # s=160 / c=140
tr=2.695  # s=1.66 / c=1.81
r_a_tr=9  # TRs for rest and activity s=8 / c=7

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
matlab -nodesktop -nosplash -r "fn_mocobatch_flex1 "${dir}

##
gnuplot "gnuplot_moco.txt"

rm ./Basis_*

################### Start here if I combine 2 runs after moco with sv_start_moco_flex.sh
################### I can combine the two runs with fslmaths, example above

# Change name to BOLD.nii
mv data_mocoBasis_a.nii BOLD.nii


### 4) Temporal filtering
echo "Temporal High-pass filter"
bptf=$(echo $block_dur/$tr | bc -l)
bptf=$(echo ${bptf%.*})
fslmaths BOLD.nii -Tmean tempmean
fslmaths BOLD.nii -bptf $bptf -1 -add tempmean BOLD_hpf.nii
# fslchfiletype NIFTI BOLD_intemp_hpf.nii.gz

#### 5) Temporal smoothing
# 3dTsmooth -prefix BOLD_hpf_sm.nii BOLD_hpf.nii

# rm BOLD_intemp.nii

#### 5) Temporal upsampling.... 
# 3dUpsample -overwrite  -datum short -prefix BOLD_intemp_hpf.nii -n 2 -input BOLD_hpf.nii
# NumVol=`3dinfo -nv BOLD_intemp.nii`
# If no just rename BOLD_hpf.nii

# cp BOLD_hpf.nii.gz BOLD_hpf.nii.gz

echo "I am correcting for the proper TR in the header"
3drefit -TR $tr BOLD_hpf.nii

## Getting the mean rest and activity, for now with LAYNII, dirty trick
# /home/laynii/LN_BOCO -Nulled BOLD_intemp_hpf.nii.gz -BOLD BOLD_intemp_hpf.nii.gz -trialBOCO $block_trs

##

echo "calculating Mean and tSNR maps"
  3dTstat  -overwrite -mean  -prefix BOLD.Mean.nii \
     BOLD_hpf.nii'[1..$]'
  3dTstat  -overwrite -cvarinv  -prefix BOLD.tSNR.nii \
     BOLD_hpf.nii'[1..$]'


# echo "curtosis and skew"
# /home/laynii/LN_SKEW -input BOLD_intemp_hpf.nii


#### 5) ACTIVTION (modify -TR_times and -stim_times in VASO and BOLD...):
# TR_times = VolumeTR(can get it from params.gen.volTR)
# stim_times = time of TR times duration of TRs in BLOCK times 2(Vaso-Bold), then just double this till number of Blocks
# stim_times = 0 TR_times*(#ON_TRs+#OFF_TRs)*2
# UBLOCK = (TR_times*#ON_TRs*2,1)
block_dur=$(echo $tr*($block_trs) | bc -l)
stim_dur=$(echo $tr*($r_a_tr) | bc -l)
# tmp=0
start=1
# end=$(echo $blocks-1 | bc -l)
# stim_times='1D: 0 '
end=$(echo $blocks-1 | bc -l)
# stim_times='1D: '
stim_times=$(echo "1D: $stim_dur ")
tmp=$stim_dur


for (( i=$start; i<=$end; i++))
do
	# tmp=$(echo $tmp+($block_dur/2) | bc -l)  # If I do temporal upsampling
	tmp=$(echo $tmp+($block_dur) | bc -l)
	stim_times=$(echo "$stim_times $tmp ")
done
tmp=$(echo $block_dur/2 | bc -l)
tmp=$(echo ${tmp%%.*})
ublock=$(echo "UBLOCK($tmp,1)")

##
echo "BOLD based on GLM"
3dDeconvolve -overwrite -jobs 16 -polort a -input BOLD_hpf.nii.gz \
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

3dTstat -mean -overwrite -prefix mean.nii BOLD_hpf.nii

3dcalc -a mean.nii -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

##### 6) Cluster activations
# Here, I need to play with the values after -1clip:
# -1clip threshold.. (~1.8), (1.5,1.2,270)
# rmm = cluster connection radius, larger value->remove small clusters
# vmul minimum cluster volume, smaller value->removes small clusters
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 3 1.4 120 2_STATS_BOLD.nii

####### 10) Masking relevant volumes
fslmaths 2_STATS_BOLD.nii -mul mask.nii clustered_BOLD_masked.nii
fslmaths BOLD.Mean.nii -mul mask.nii BOLD_mean_masked.nii

### Create the manual Visual and Motor mask....
##### 11) Masking again for a ROI, if needed....
fslmaths mask.nii -mul roi_visual_motor.nii.gz roi_visual_motor_new.nii
fslmaths clustered_BOLD_masked.nii -mul roi_visual_motor_new.nii clustered_BOLD_masked.nii


##### 7) Get mean of activations and percentage change
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
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode Cu -overwrite -prefix scaled_BOLD.nii -input clustered_BOLD_masked.nii.gz 


### NOW CREATE Manually the mask in FSLEYES, call it rim.nii

#estimating layers based on rim
LN_GROW_LAYERS -rim rim.nii -N 11 -vinc 40

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

# BOLD
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,3,5)' -prefix BOLD_3_5.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,5.001,7)' -prefix BOLD_5_7.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,7.001,9)' -prefix BOLD_7_9.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,9.001,11)' -prefix BOLD_9_11.nii
3dcalc -overwrite -a clustered_BOLD_masked.nii -expr 'within(a,11.001,50)' -prefix BOLD_11_plus.nii


# Calculate volume mean tSNR
3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./BOLD_tSNR.nii






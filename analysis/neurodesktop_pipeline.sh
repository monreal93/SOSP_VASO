# Neurodesktop tools to use:
ml afni
ml fsl
ml laynii

cd /neurodesktop-storage/5T3/Alejandro/sosp_vaso/data

folder="08282024_sb_9T"
scan="cb_01"
r_a_tr=15     # rest/activity TRs paper=sv(8), cv(7)
tr=1.69  # volume TR paper=sv(1.66), cv(1.81)
recon_file="cb_01_b_epi"  # name of the reconstructred volume.. it should be in the recon folder...
motion_glm=1     # GLM including motion parametrs

cd ${folder}

mkdir ./analysis/${scan}
chmod ugo+rwx ./analysis/${scan}
cd ./analysis/${scan}

# gre1=../../tmp/${scan}.nii
file=../../recon/${recon_file}.nii

vol=$(3dinfo -nv ${file})

blocks=$(echo $vol/$r_a_tr/2 | bc -l)
blocks=$(echo ${blocks%.*})
block_dur=$(echo $tr*$blocks | bc -l)
block_dur=$(echo ${block_dur%.*})
block_trs=$(echo $r_a_tr*2 | bc -l)
block_trs=$(echo ${block_trs%.*})

# ToDo: Find an easy way to realign Anatomy as fMRI...
# Center and deoblique datasets
3dinfo -obliquity ${file} >> obliquity.txt 
3drefit -xorigin cen -yorigin cen -zorigin cen ${file}
# 3drefit -xorigin cen -yorigin cen -zorigin cen ${gre1}
3dWarp -deoblique ${file}
# 3dWarp -deoblique ${gre1}

# ) Making sure orientation is LPI
3drefit -orient LPI ${file}
# 3drefit -orient LPI ${gre1}

# # # Sometimes I need to Flip in y direction
# 3dLRflip -Y -prefix ${file} ${file} -overwrite
# 3dLRflip -Y -prefix ${gre1} ${gre1} -overwrite

# 1) Creating Mask
#--- with FSL
# bet $gre1 ./tmp -Z -f 0.4 -g 0 -n -m
# mv tmp_mask.nii.gz mask.nii.gz
# gzip -d mask.nii.gz
#--- with AFNI

# if [ "${scan:0:1}" = "s" ]; then
#     3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${gre1} -overwrite
# else
    3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${file} -overwrite
# fi

3drefit -xdel $(3dinfo -adi ${file}) mask.nii
3drefit -ydel $(3dinfo -adj ${file}) mask.nii
3drefit -zdel $(3dinfo -adk ${file}) mask.nii

# Replacing non-steady state volumes..
# ToDo: How many volumes should I replace?
3dTcat -prefix ${file} ${file}'[4..7]' ${file}'[4..$]' -overwrite

# ) Motion correction
# ToDo: Check if mask can be used...
# ToDo: Use same function as Renzo's new scripts... 3dAllineate
echo "Motion correction..."
# --- 3dvolreg
3dvolreg -base 4 -heptic -zpad 1 -overwrite -prefix ./${scan}_mc.nii -1Dfile motion.1D ${file}
# --- 3dAllineate
# 3dTstat -mean -overwrite -prefix ./n_ref.nii ${b_file}'[0..3]'  # create reference
# 3dAllineate -1Dmatrix_save  matrix.aff12.1D -1Dparam_save   param.aff12 -cost lpa \
#     -prefix ./${scan}_b_mc.nii -base ${output_dir} ./n_ref.nii -source ${b_file} \
#     -weight ./mask.nii -warp shift_rotate -final wsinc5

3drefit -TR $tr ./${scan}_mc.nii 

#### ) Temporal filtering
# This commands can be used for Tmporal filtering in AFNI 3dDetrend, 3dBandpass
echo "Temporal High-pass filter..."
bptf=$(echo $block_dur/$tr | bc -l)
bptf=$(echo ${bptf%.*})
fslmaths ./${scan}_mc.nii -Tmean tempmean
fslmaths ./${scan}_mc.nii -bptf $bptf -1 -add tempmean ./${scan}_mc_hpf.nii
fslchfiletype NIFTI ${scan}_mc_hpf.nii

#### ) Quality metrics
echo "calculating Mean and tSNR maps ..."
3dTstat -mean -prefix mean.nii ./${scan}_mc_hpf.nii'[1..$]' -overwrite
3dTstat  -overwrite -cvarinv  -prefix tSNR.nii ./${scan}_mc_hpf.nii'[1..$]'

#### ) Activation maps
block_dur=$(echo $tr*$block_trs | bc -l)
block_dur=$(echo ${block_dur%.*})
tmp=0
start=1
end=$(echo $blocks-1 | bc -l) # original
tmp=$(echo $block_dur/2 | bc -l)
tmp=$(echo ${tmp%%.*})
ublock=$(echo "UBLOCK($tmp,1)")
# stim_times='1D: 0 '  # original 
# end=$(echo $blocks-2 | bc -l) # new
stim_times=$(echo "1D: $tmp ") # new
# tmp=$(echo $block_dur | bc -l)  # If temporal upsampling # New
# stim_times=$(echo "$stim_times $tmp ") # new

for (( i=$start; i<=$end; i++))
do
    tmp=$(echo $tmp+$block_dur | bc -l)
	stim_times=$(echo "$stim_times $tmp ")
done

if [ "$motion_glm" = 0 ]; then
    #### BOLD
    echo "BOLD based on GLM..."
    3dDeconvolve -overwrite -jobs 16 -polort a -input ./${scan}_mc_hpf.nii\
                -num_stimts 1 \
                -TR_times $tr \
                -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
                -tout \
                -x1D MODEL_wm \
                -iresp 1 HRF_BOLD.nii \
                -errts residual_BOLD.nii \
                -bucket STATS_BOLD.nii
else
    #### BOLD with Motion....
    echo "BOLD based on GLM..."
    3dDeconvolve -overwrite -jobs 16 -polort a \
                -force_TR $tr \
                -input ./${scan}_mc.nii\
                -num_stimts 7 \
                -TR_times $tr \
                -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
                -tout \
                -x1D MODEL_wm \
                -iresp 1 HRF_BOLD.nii \
                -errts residual_BOLD.nii \
                -stim_file 2 motion.1D'[0]' -stim_base 2 -stim_label 2 roll \
                -stim_file 3 motion.1D'[1]' -stim_base 3 -stim_label 3 pitch \
                -stim_file 4 motion.1D'[2]' -stim_base 4 -stim_label 4 yaw \
                -stim_file 5 motion.1D'[3]' -stim_base 5 -stim_label 5 dS \
                -stim_file 6 motion.1D'[4]' -stim_base 6 -stim_label 6 dL \
                -stim_file 7 motion.1D'[5]' -stim_base 7 -stim_label 7 dP \
                -bucket STATS_BOLD.nii
fi

3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a'    -prefix 1_HRF_BOLD.nii   -overwrite 
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_BOLD.nii   -overwrite 

3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'    -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_BOLD.nii -overwrite
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr 'a'    -prefix 2_STATS_BOLD.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii ./${scan}_mc_hpf.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

### BOLD mean residual
3dcalc -overwrite -a clustered_BOLD.nii -expr 'step(a-1.8)' -prefix b_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask b_bin_output.nii -1DRformat -quiet ./residual_BOLD.nii > residual_BOLD.dat

####### ) Masking relevant volumes
3dcalc -a mean.nii -b mask.nii -expr 'a*b' -prefix mean_msk.nii -overwrite
3dcalc -a 2_STATS_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii -overwrite
3dcalc -a tSNR.nii -b mask.nii -expr 'a*b' -prefix tSNR_msk.nii -overwrite

# Let's ommit the first two slices (fold over artifacts)
# ToDo: Find better way of doing this
slices=$(3dinfo -nk BOLD_msk.nii)
slices=$(($slices-1))
# BOLD
3dZcutup -keep 2 "${slices}" -prefix BOLD_msk.nii BOLD_msk.nii -overwrite
3dZeropad -I 2 -prefix BOLD_msk.nii BOLD_msk.nii -overwrite

##### ) Cluster activations
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 3 1.4 120 BOLD_msk.nii

##### ) Get mean of activations and percentage change
### BOLD
3dcalc -overwrite -a clustered_BOLD.nii -expr 'step(a-1.8)' -prefix bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask bin_output.nii -1DRformat -quiet ./${scan}_mc_hpf.nii > timecourse_BOLD.dat

####### ) Mean tSNR and effective tSNR
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./tSNR_msk.nii) 
echo -e "BOLD brain mean tSNR: \n $mean_tSNR" >> results.txt

3dcalc -a ./tSNR_msk.nii -expr "a/sqrt($tr)" -prefix ./eff_tSNR.nii -overwrite
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./eff_tSNR.nii) 
echo -e "BOLD brain mean effective tSNR: \n $mean_tSNR_b" >> results.txt



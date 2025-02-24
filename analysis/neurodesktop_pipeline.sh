# Neurodesktop tools to use:
ml afni
ml fsl
ml laynii

cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="02202025_sb_9T_paper"
scan="sb_601_DS_SO_06mm_18fovz_12te_6te"
traj="girf"   # nom or girf
r_a_tr=7   # rest/activity TRs
tr=2.7 # volume TR paper=sv(1.66), cv(1.81)
recon_file="sb_601_DS_SO_06mm_18fovz_12te_6te_ech2_girf_cs_fsb0_co_rDORK"  # name of the reconstructred volume.. it should be in the recon folder...
motion_glm=0     # GLM including motion parametrs
stim=1           # number of stimulus for block desig

cd ${folder}

mkdir ./analysis/${scan}_${traj}
chmod ugo+rwx ./analysis/${scan}_${traj}
cd ./analysis/${scan}_${traj}

### Starting... 
gre1=../../tmp/${scan}_${traj}.nii
file=../../recon/${recon_file}.nii
file_rev_ph=../../recon/${recon_file}_rev_ph.nii  # Reversed phase enc EPI volume...

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
3drefit -xorigin cen -yorigin cen -zorigin cen ${gre1}
3dWarp -deoblique ${file}
3dWarp -deoblique ${gre1}

# ) Making sure orientation is LPI (Transversal) or RIA (Coronal)
3drefit -orient RIA ${file}
3drefit -orient RIA ${gre1}

# # # Sometimes I need to Flip in y direction
# 3dLRflip -Y -prefix ${file} ${file} -overwrite
# 3dLRflip -Y -prefix ${gre1} ${gre1} -overwrite

# If EPI, then do TOPUP distorsion correction (WIP)
# if [ "${scan:0:1}" = "c" ]; then
    # 3drefit -xorigin cen -yorigin cen -zorigin cen ${file_rev_ph}
    # 3dWarp -deoblique ${file_rev_ph}
    # 3drefit -orient LPI ${file_rev_ph}
    # 3dTcat -prefix ${scan}_topup.nii ${file_rev_ph}'[2]' ${file} -overwrite # concatenate 2nd volume of reversed PE and fMRI series

    # # AMM: ToDo: Find good way to crete this text file.. 
    # printf "0    -1    0    0.01" > acq_params.txt
    # for ((i = 0 ; i < $vol ; i++ ));
    # do
    #     printf "\n0    1    0    0.01" >> acq_params.txt
    # done
    # echo "Topup correction..."
    # topup --imain=${scan}_topup.nii --datain=acq_params.txt --config=b02b0.cnf --iout=${file}_topup1.nii --nthr=8 --fwhm=0

    # 3dTcat -prefix ${scan}_topup.nii ${scan}_topup1.nii'[1..$]' -overwrite # concatenate 2nd volume of reversed PE and fMRI series
    # rm ${scan}_topup1.nii
    # file=${scan}_topup.nii
# fi

# 1) Creating Mask
#--- with FSL
bet $gre1 ./tmp -Z -f 0.5 -g 0 -n -m
mv tmp_mask.nii.gz mask.nii.gz
gzip -d mask.nii.gz
3drefit -xorigin cen -yorigin cen -zorigin cen mask.nii
#--- with AFNI
# if [ "${scan:0:1}" = "s" ]; then
#     3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${gre1}'[2]' -overwrite
# else
#     3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${file} -overwrite
# fi

# 3drefit -xdel $(3dinfo -adi ${file}) mask.nii
# 3drefit -ydel $(3dinfo -adj ${file}) mask.nii
# 3drefit -zdel $(3dinfo -adk ${file}) mask.nii

# Replacing non-steady state volumes..
# ToDo: How many volumes should I replace?
3dTcat -prefix ${file} ${file}'[4..7]' ${file}'[4..$]' -overwrite

# ) Motion correction
# ToDo: Check if mask can be used...
# ToDo: Use same function as Renzo's new scripts... 3dAllineate
echo "Motion correction..."
# --- 3dvolreg, only works for small motion...
# 3dvolreg -base 4 -cubic -zpad 1 -overwrite -prefix ./${scan}_mc.nii -1Dfile motion.1D ${file}
# --- 3dAllineate
3dTstat -mean -overwrite -prefix ./n_ref.nii ${file}'[2..3]'  # create reference
3dAllineate -1Dmatrix_save  matrix.aff12.1D -1Dparam_save   param.aff12 -cost lpa \
    -prefix ./${scan}_mc.nii -base ${output_dir} ./n_ref.nii -source ${file} \
    -weight ./mask.nii -warp shift_rotate -final wsinc5 -overwrite

3drefit -TR $tr ./${scan}_mc.nii 

#### ) Temporal filtering
# This commands can be used for Tmporal filtering in AFNI 3dDetrend, 3dBandpass
echo "Temporal High-pass filter..."
bptf=$(echo $block_dur/$tr | bc -l)
bptf=$(echo ${bptf%.*})
fslmaths ./${scan}_mc.nii -Tmean tempmean.nii.gz
fslmaths ./${scan}_mc.nii -bptf $bptf -1 -add tempmean.nii.gz ./${scan}_mc_hpf.nii
fslchfiletype NIFTI ${scan}_mc_hpf.nii

#### ) Quality metrics
echo "calculating Mean and tSNR maps ..."
3dTstat -mean -prefix mean.nii ./${scan}_mc_hpf.nii'[1..$]' -overwrite
3dTstat -cvarinv  -prefix tSNR.nii ./${scan}_mc_hpf.nii'[1..$]' -overwrite

#### ) Activation maps
block_dur=$(echo $tr*$block_trs | bc -l)
block_dur=$(echo ${block_dur%.*})
tmp=0
start=1
dummy=0
end=$(echo $blocks-1 | bc -l) # original
tmp=$(echo $block_dur/2 | bc -l)
tmp=$(echo ${tmp%%.*})
ublock=$(echo "UBLOCK($tmp,1)")

stim_times=$(echo "1D: $tmp ") # new

if [ "$stim" = 2 ]; then
    tmp1=$(echo $block_dur+$block_dur/2 | bc -l)
    tmp1=$(echo ${tmp1%%.*})
    stim_times1=$(echo "1D: $tmp1 ") # new
    dummy=$block_dur
fi

# Adding timing for all blocks...
for (( i=$start; i<=$end; i++))
do
    tmp=$(echo $tmp+$block_dur+$dummy | bc -l)
    tmp1=$(echo $tmp1+$block_dur+$dummy | bc -l)
	stim_times=$(echo "$stim_times $tmp ")
    stim_times1=$(echo "$stim_times1 $tmp1 ")
done

if [ "$motion_glm" = 0 ]; then
    #### BOLD
    echo "BOLD based on GLM..."
    if [ "$stim" = 2 ]; then
        3dDeconvolve -overwrite -jobs 16 -polort a -input ./${scan}_mc_hpf.nii\
                    -num_stimts 2 \
                    -TR_times $tr \
                    -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
                    -stim_times 2 "$stim_times1" "$ublock" -stim_label 2 Task1 \
                    -tout \
                    -x1D MODEL_wm \
                    -iresp 1 HRF_BOLD.nii \
                    -errts residual_BOLD.nii \
                    -bucket STATS_BOLD.nii
    else
        3dDeconvolve -overwrite -jobs 16 -polort a -input ./${scan}_mc_hpf.nii\
            -num_stimts 1 \
            -TR_times $tr \
            -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
            -tout \
            -x1D MODEL_wm \
            -iresp 1 HRF_BOLD.nii \
            -errts residual_BOLD.nii \
            -bucket STATS_BOLD.nii
    fi
else
    #### BOLD with Motion....
    echo "BOLD based on GLM..."
    if [ "$stim" = 2 ]; then
        3dDeconvolve -overwrite -jobs 16 -polort a \
                    -force_TR $tr \
                    -input ./${scan}_mc_hpf.nii\
                    -num_stimts 8 \
                    -TR_times $tr \
                    -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
                    -stim_times 2 "$stim_times1" "$ublock" -stim_label 2 Task1 \
                    -tout \
                    -x1D MODEL_wm \
                    -iresp 1 HRF_BOLD.nii \
                    -errts residual_BOLD.nii \
                    -stim_file 3 motion.1D'[0]' -stim_base 3 -stim_label 3 roll \
                    -stim_file 4 motion.1D'[1]' -stim_base 4 -stim_label 4 pitch \
                    -stim_file 5 motion.1D'[2]' -stim_base 5 -stim_label 5 yaw \
                    -stim_file 6 motion.1D'[3]' -stim_base 6 -stim_label 6 dS \
                    -stim_file 7 motion.1D'[4]' -stim_base 7 -stim_label 7 dL \
                    -stim_file 8 motion.1D'[5]' -stim_base 8 -stim_label 8 dP \
                    -bucket STATS_BOLD.nii
    else
        3dDeconvolve -overwrite -jobs 16 -polort a \
            -force_TR $tr \
            -input ./${scan}_mc_hpf.nii\
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
fi

3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a'    -prefix 1_HRF_BOLD.nii   -overwrite 
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_BOLD.nii   -overwrite 

3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'    -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_BOLD.nii -overwrite
# Task 1
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr 'a'    -prefix 2_STATS_BOLD.nii -overwrite 

# Task 2
3dcalc -a STATS_BOLD.nii'[4]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD_2.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[4]'  -expr 'a'    -prefix 2_STATS_BOLD_2.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii ./${scan}_mc_hpf.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii


####### ) Masking relevant volumes
3dcalc -a mean.nii -b mask.nii -expr 'a*b' -prefix mean_msk.nii -overwrite
3dcalc -a tSNR.nii -b mask.nii -expr 'a*b' -prefix tSNR_msk.nii -overwrite
3dcalc -a 2_STATS_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii -overwrite

# Let's ommit the first two slices (fold over artifacts)
# ToDo: Find better way of doing this
slices=$(3dinfo -nk BOLD_msk.nii)
slices=$(($slices-1))
# BOLD
3dZcutup -keep 2 "${slices}" -prefix BOLD_msk.nii BOLD_msk.nii -overwrite
3dZeropad -master mask.nii -prefix BOLD_msk.nii BOLD_msk.nii -overwrite

##### ) Cluster activations, sometimes I use -sided RIGHT_TAIL 3 , sometimes 4 (check the clustering result...)
# 3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 3 1.4 120 BOLD_msk.nii
thr=1.6
3dClusterize -inset BOLD_msk.nii -ithr 0 -idat 0 -NN 1 -1sided RIGHT_TAIL $thr -clust_nvox 40 -pref_dat clustered_BOLD.nii -overwrite

########### Writing some metrics in results.txt
### Percentage of active voxels
active_vox=$(3dBrickStat -count -non-zero -mask ./mask.nii ./clustered_BOLD.nii)
mask_vox=$(3dBrickStat -count -non-zero -mask ./mask.nii ./mask.nii)
active_vox_per=$(echo $active_vox/$mask_vox*100 | bc -l)
echo -e "Percentage active voxels threshold = $thr \n $active_vox_per" >> results.txt

### Median z-score
med_z_score=$(3dBrickStat -median -non-zero -mask ./mask.nii ./clustered_BOLD.nii)
echo -e "Median z-score: \n $med_z_score" >> results.txt

####### ) Mean tSNR and effective tSNR
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./tSNR_msk.nii) 
echo -e "BOLD brain mean tSNR: \n $mean_tSNR_b" >> results.txt

## Effective tSNR
3dcalc -a ./tSNR_msk.nii -expr "a/sqrt($tr)" -prefix ./eff_tSNR.nii -overwrite
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./eff_tSNR.nii) 
echo -e "BOLD brain mean effective tSNR: \n $mean_tSNR_b" >> results.txt

###########

# ## temporal Contrast to Noise Ratio
# 3dcalc -prefix b_bin_output_neg.nii -a mask.nii -b b_bin_output.nii -expr "a-b" -overwrite
# std_activity=$(3dBrickStat -stdev -non-zero -mask b_bin_output.nii ./${scan}_mc_hpf.nii)
# std_noise=$(3dBrickStat -stdev -non-zero -mask b_bin_output_neg.nii ./${scan}_mc_hpf.nii)
# # 3dTstat -prefix std_act.nii -mask b_bin_output.nii -stdev ./${scan}_mc_hpf.nii -overwrite -overwrite
# # 3dTstat -prefix std_noise.nii -mask b_bin_output_neg.nii -stdev ./${scan}_mc_hpf.nii -overwrite -overwrite
# # std_activity=$(3dBrickStat -stdev -non-zero -mask b_bin_output.nii std_act.nii)
# # std_noise=$(3dBrickStat -stdev -non-zero -mask b_bin_output_neg.nii std_noise.nii)
# tCNR_b=$(echo $std_activity/$std_noise | bc -l)
# echo -e "tCNR: \n $tCNR_b" >> results.txt

##### ) Get mean of activations and percentage change
### BOLD
3dcalc -overwrite -a clustered_BOLD.nii -expr "step(a-$thr)" -prefix bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask bin_output.nii -1DRformat -quiet ./${scan}_mc_hpf.nii > timecourse_BOLD.dat

### BOLD mean residual
3dcalc -overwrite -a clustered_BOLD.nii -expr "step(a-$thr)" -prefix b_bin_output.nii
# Might be needed to clean up the binary mask before getting mean
3dROIstats -mask b_bin_output.nii -1DRformat -quiet ./residual_BOLD.nii > residual_BOLD.dat

################################## ROI quality metrics....
# ToDo: Open clustered_BOLD and mean_mask.nii in FSLeyes, then create a mask with brush size 6 covering the desired area...
### Percentage of active voxels
# thr=1.6
active_vox=$(3dBrickStat -count -non-zero -mask ./roi_mask.nii ./clustered_BOLD.nii)
mask_vox=$(3dBrickStat -count -non-zero -mask ./roi_mask.nii ./roi_mask.nii)
active_vox_per=$(echo $active_vox/$mask_vox*100 | bc -l)
echo -e "ROI Percentage active voxels threshold = $thr \n $active_vox_per" >> results.txt

### Median z-score
3dcalc -overwrite -a clustered_BOLD.nii -b roi_mask.nii -expr "a*b" -prefix clustered_BOLD_ROI.nii
med_z_score=$(3dBrickStat -median -non-zero -mask ./roi_mask.nii ./clustered_BOLD.nii)
echo -e "ROI Median z-score: \n $med_z_score" >> results.txt

####### ) Mean tSNR and effective tSNR
3dcalc -overwrite -a tSNR_msk.nii -b roi_mask.nii -expr "a*b" -prefix tSNR_ROI.nii
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./roi_mask.nii ./tSNR_msk.nii) 
echo -e "ROI BOLD brain mean tSNR: \n $mean_tSNR_b" >> results.txt
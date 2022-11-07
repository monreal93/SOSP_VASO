
# Make sure you change the directory to the one where you are running the analysis....
# sosp_vaso/data/"date_of_scan"/analysis/...

# If I want to combine 2 diff runs, I can do it with fslmaths, after moco I can average them, example:
# fslmaths ./Documents/PhD/PhD_2022/sosp_vaso/data/_10282022_sv/analysis/sv_02/Nulled_Basis_b.nii -add ./Documents/PhD/PhD_2022/sosp_vaso/data/_10282022_sv/analysis/sv_03/Nulled_Basis_b.nii -div 2 ./Documents/PhD/PhD_2022/sosp_vaso/data/_10282022_sv/analysis/sv_02_03/Nulled_Basis_b.nii 

#### 1) Let's first copy the files I need into this folder
cp /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/analysis/youtube/gnuplot_moco.txt ./gnuplot_moco.txt

# I am not sure if I also want to copy mocobatch.m
# cp /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/analysis/youtube/

#### 2) MOTION CORRECTION (modify cp path in mocobatch.m):
path=$(ls -1 ../cv_*.nii | head -n 1)
dir=$(pwd)

3dAutomask -prefix moma.nii -peels 3 -dilate 2 ${path}
#3dAutomask -prefix moma.nii -peels 3 -dilate 2 *v_0*.nii 

cp ./*v_0* ./Basis_a.nii
3dTcat -prefix Basis_a.nii Basis_a.nii'[4..7]' Basis_a.nii'[4..$]' -overwrite
cp ./Basis_a.nii ./Basis_b.nii

3dinfo -nt Basis_a.nii >> NT.txt
3dinfo -nt Basis_b.nii >> NT.txt

###### 3) STOP HEREE...... OPEN MATLAB..... ############
###### For now, open matlab and run file mocobatch.m and make sure you change the path cd... in the .m file.... ########
###### STOP HEREE...... OPEN MATLAB..... ############
matlab -nodesktop -nosplash -r "fn_mocobatch_flex "${dir}

gnuplot "gnuplot_moco.txt"

# rm ./Basis_*

################### Start here if I combine 2 runs after moco with sv_start_moco_flex.sh
################### I can combine the two runs with fslmaths, example above

#### 4) BOLD CORRECTION (modify -trialBOCO, and 3drefit -TR):
3dcalc -a Nulled_Basis_b.nii'[1..$(2)]' -expr 'a' -prefix Nulled.nii -overwrite
3dcalc -a Not_Nulled_Basis_a.nii'[0..$(2)]' -expr 'a' -prefix BOLD.nii -overwrite
3dUpsample -overwrite  -datum short -prefix Nulled_intemp.nii -n 2 -input Nulled.nii
3dUpsample -overwrite  -datum short -prefix BOLD_intemp.nii   -n 2 -input   BOLD.nii
NumVol=`3dinfo -nv BOLD_intemp.nii`
3dTcat -overwrite -prefix Nulled_intemp.nii Nulled_intemp.nii'[0]' Nulled_intemp.nii'[0..'`expr $NumVol - 2`']' 

/home/laynii/LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO 18

## 
echo "I am correcting for the proper TR in the header"
3drefit -TR 1.67 BOLD_intemp.nii
3drefit -TR 1.67 VASO_LN.nii

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

echo "VASO based on GLM" 
3dDeconvolve -overwrite -jobs 16 -polort a -input VASO_LN.nii\
             -num_stimts 1 \
             -TR_times 1.67 \
             -stim_times 1 '1D: 0 53.4 106.8 160.2 267 320.4 373.8 427.2 480.6' 'UBLOCK(26.7,1)' -stim_label 1 Task \
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

##
echo "BOLD based on GLM"
3dDeconvolve -overwrite -jobs 16 -polort a -input BOLD_intemp.nii\
             -num_stimts 1 \
             -TR_times 1.67 \
             -stim_times 1 '1D: 0 53.4 106.8 160.2 267 320.4 373.8 427.2 480.6' 'UBLOCK(26.7,1)' -stim_label 1 Task \
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
3dcalc -a mean.nii -b  1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

##### 6) Cluster activations
3dclust -1noneg -overwrite -prefix clustered_VASO.nii -1clip 1.5 1.2 270 2_STATS_VASO.nii

##### 7) Get mean of activations
3dcalc -overwrite -a activationmap.nii -expr 'step(a-1.8)' -prefix binarised_output.nii

3dROIstats -mask binarised_output.nii -1DRformat -quiet time_series.nii > timecourse.dat



#### 7) Possible smoothing
LN_GRADSMOOTH -input 2_STATS_VASO.nii -gradfile T1_weighted.nii -FWHM 1.5 -within -selectivity 0.08 







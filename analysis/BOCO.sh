#!/bin/bash


echo "The script starts now:  I expect two files Not_Nulled_Basis_a.nii and Nulled_Basis_b.nii that are motion corrected with SPM"


echo "temporal upsampling and shifting happens now"
3dcalc -a Nulled_Basis_b.nii'[1..$(2)]' -expr 'a' -prefix Nulled.nii -overwrite
3dcalc -a Not_Nulled_Basis_a.nii'[0..$(2)]' -expr 'a' -prefix BOLD.nii -overwrite
3dUpsample -overwrite  -datum short -prefix Nulled_intemp.nii -n 2 -input Nulled.nii
3dUpsample -overwrite  -datum short -prefix BOLD_intemp.nii   -n 2 -input   BOLD.nii
NumVol=`3dinfo -nv BOLD_intemp.nii`
3dTcat -overwrite -prefix Nulled_intemp.nii Nulled_intemp.nii'[0]' Nulled_intemp.nii'[0..'`expr $NumVol - 2`']' 

echo "BOLD correction happens now"
#/home/laynii/LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO 40
/home/laynii/LN_BOCO -Nulled Nulled_intemp.nii -BOLD BOLD_intemp.nii -trialBOCO 10

echo "I am correcting for the proper TR in the header"
# 3drefit -TR 1.5 BOLD_intemp.nii
# 3drefit -TR 1.5 VASO_LN.nii
3drefit -TR 1.3589 BOLD_intemp.nii
3drefit -TR 1.3589 VASO_LN.nii

echo "calculating T1 in EPI space"
NumVol=`3dinfo -nv Nulled_Basis_b.nii`
3dcalc -a Nulled_Basis_b.nii'[3..'`expr $NumVol - 2`']' -b  Not_Nulled_Basis_a.nii'[3..'`expr $NumVol - 2`']' -expr 'a+b' -prefix combined.nii -overwrite
3dTstat -cvarinv -prefix T1_weighted.nii -overwrite combined.nii 
rm combined.nii

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





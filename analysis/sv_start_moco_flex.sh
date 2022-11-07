#!/bin/bash

cnt=0

path=$(ls -1 ../sv_*.nii | head -n 1)
dir=$(pwd)

3dAutomask -prefix moma.nii -peels 3 -dilate 2 ${path} 

# 3dAutomask -prefix moma.nii -peels 3 -dilate 2 S*.nii 

for filename in ../sv_0*.nii
do
echo $filename
cp $filename ./Basis_${cnt}a.nii
3dTcat -prefix Basis_${cnt}a.nii Basis_${cnt}a.nii'[4..7]' Basis_${cnt}a.nii'[4..$]' -overwrite
cp ./Basis_${cnt}a.nii ./Basis_${cnt}b.nii

3dinfo -nt Basis_${cnt}a.nii >> NT.txt
3dinfo -nt Basis_${cnt}b.nii >> NT.txt
cnt=$(($cnt+1))

done

###### STOP HEREE...... OPEN MATLAB..... ############
###### For now, open matlab and run file mocobatch.m and make sure you change the path cd... in the .m file.... ########
###### STOP HEREE...... OPEN MATLAB..... ############

matlab -nodesktop -nosplash -r "fn_mocobatch_flex "${dir}

#### Let's first copy the files I need into this folder
cp /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/analysis/Alejandro/gnuplot_moco2.txt ./gnuplot_moco2.txt
gnuplot "gnuplot_moco2.txt"

rm ./Basis_*.nii

# Now lets get the mean of both runs
3dMean -prefix Nulled_Basis_b.nii Nulled_Basis_0b.nii Nulled_Basis_1b.nii
3dMean -prefix Not_Nulled_Basis_a.nii Not_Nulled_Basis_0a.nii Not_Nulled_Basis_1a.nii

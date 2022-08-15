#!/bin/bash

cd /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/data/08052022_sv_abc/analysis/sv_01

# cp Smagn.nii ./Basis_a.nii
cp ./sv_01_v_b-b_v.nii ./Basis_a.nii
3dTcat -prefix Basis_a.nii Basis_a.nii'[4..7]' Basis_a.nii'[4..$]' -overwrite
cp ./Basis_a.nii ./Basis_b.nii

3dinfo -nt Basis_a.nii >> NT.txt
3dinfo -nt Basis_b.nii >> NT.txt


# NOTE: matlab might be at a different loation on your mashine. 
# /Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nosplash -r "mocobatch"

matlab -nodesktop -nosplash 

# -r "/home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/analysis/renzo_example/EPI_analysis/mocobatch.m"

# AMM: In MATLAB RUN
# run('./sosp_vaso/analysis/renzo_example/EPI_analysis/mocobatch.m')

gnuplot "gnuplot_moco.txt"

rm ./Basis_*.nii

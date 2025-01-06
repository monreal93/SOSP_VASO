# After upsampling, you need to create a rim.nii mask..
# follow video: https://www.youtube.com/watch?v=VUDi2Iskzz4&t=913s
# starting from 26:00
# Neurodesktop tools to use:
ml afni
ml laynii

cd /neurodesktop-storage/5T3/Alejandro/sosp_vaso/data

folder="08302024_sb_9T"
scan="sb_22_OUT_2shot_6te_14fa"
t1_file="sb_22_OUT_2shot_6te_14fa_t1_full_brain"

cd ${folder}$"/analysis/"${scan}

# if [ "${scan:1:2}" = "v" ]; then
#     t1_file=T1_msk.nii # FOR VASO
# else
#     t1_file=${scan}_t1.nii # FOR BOLD ...
# fi

##### 9) Layering
# Upsample
delta_x=$(3dinfo -di ${t1_file}.nii)
delta_y=$(3dinfo -dj ${t1_file}.nii)
delta_z=$(3dinfo -dk ${t1_file}.nii)
sdelta_x=$(echo "((sqrt($delta_x * $delta_x) / 5))"|bc -l)
sdelta_y=$(echo "((sqrt($delta_y * $delta_y) / 5))"|bc -l)
sdelta_z=$(echo "((sqrt($delta_z * $delta_z) / 1))"|bc -l) 

## Scaling with AFNI... Not the best for oblique datasets
# here I only upscale in 2 dimensions. 
#3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode NN -overwrite -prefix scaled_$1 -input $1 
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode Cu -overwrite -prefix scaled_T1.nii -input ${t1_file}.nii
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode Cu -overwrite -prefix scaled_VASO.nii -input clustered_VASO.nii
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z -rmode Cu -overwrite -prefix scaled_BOLD.nii -input clustered_BOLD.nii


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

# 1dplot -sepscl layer.dat 
mv layer.dat layer_VASO.dat

#extractiong profiles for BOLD
#get mean value, STDEV, and number of voxels
3dROIstats -mask rim_layers.nii -1DRformat -quiet -nzmean scaled_BOLD.nii > layer_t.dat
3dROIstats -mask rim_layers.nii -1DRformat -quiet -sigma scaled_BOLD.nii >> layer_t.dat
3dROIstats -mask rim_layers.nii -1DRformat -quiet -nzvoxels scaled_BOLD.nii >> layer_t.dat
#format file to be in columns, so gnuplot can read it.
WRD=$(head -n 1 layer_t.dat|wc -w); for((i=2;i<=$WRD;i=i+2)); do awk '{print $'$i'}' layer_t.dat| tr '\n' ' ';echo; done > layer.dat

# 1dplot -sepscl layer.dat 
mv layer.dat layer_BOLD.dat
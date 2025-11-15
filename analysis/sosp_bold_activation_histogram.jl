using Pkg

cd("/neurodesktop-storage/5T4/Alejandro/sosp_vaso/")
Pkg.activate("./analysis/")

using NIfTI
using CairoMakie
using MIRTjim: jim, prompt; jim(:prompt, true)
using Infiltrator
using Colors
using Statistics


directories = [
                "01152025_sb_9T_paper",
                "01152025_sb_9T_paper",
                "02062025_sb_9T_paper",
                "02062025_sb_9T_paper",
                "02132025_sb_9T_paper",
                "02132025_sb_9T_paper",
                "02202025_sb_9T_paper",
                ]

# List all directories.. it will plot 1-2 per each two-sided violin plot
scans = [
        "sb_001_DS_SO_06mm_18fovz_12te_6te",
        "sb_001_DS_SO_06mm_18fovz_12te_6te",
        "sb_011_DS_SO_06mm_18fovz_12te_6te",
        "sb_011_DS_SO_06mm_18fovz_12te_6te",
        "sb_311_DS_SO_06mm_18fovz_12te_6te",
        "sb_311_DS_SO_06mm_18fovz_12te_6te",
        "sb_312_DS_SO_06mm_18fovz_12te_6te",
        "sb_312_DS_SO_06mm_18fovz_12te_6te",
        "sb_121_DS_SO_06mm_18fovz_12te_6te",
        "sb_121_DS_SO_06mm_18fovz_12te_6te",
        "sb_122_DS_SO_06mm_18fovz_12te_6te",
        "sb_122_DS_SO_06mm_18fovz_12te_6te",
        "sb_601_DS_SO_06mm_18fovz_12te_6te",
        "sb_601_DS_SO_06mm_18fovz_12te_6te",
        ]


drive = ["5T4","5T4","5T4","5T4","5T4","5T4","5T4"]

# echo=["te2","te1"] 
echo=["ech2","ech1"] 

# colors=[RGBA{Float32}(1,0,0,0.3),RGBA{Float32}(0,0,1,0.3)]
colors = ["(:red,0.5)", "(:blue,0.5)"]

plot_range=[100,100]  # [30,50] for tSNR, []

#####
# tSNR distribution in GM over subjects (Violin plots)
####### BOLD/VASO
fig = Figure(size = (1200,500))         # For 3 "subjects"
# ax = Axis(fig[1, 1])
ax = Axis(fig[1,1], limits = (0.5,18,-0.1,1.8), yticklabelsize=64)
data = Vector{Float32}(undef,0)
category =  Vector{Float32}(undef,0)
# color = Vector{RGBA{Float32}}(undef,0)
dodge = Vector{Int64}(undef,0)
side = Vector{Symbol}(undef,0)
@infiltrate

for kk=1:1
    p=0
    for i=1:(length(directories))
        # for p=1:length(scans)
                for k=1:length(echo)
                    p +=1
                    #only common active voxels
                    # tmp_data=niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p][1:6],"_te_comparison","/",scans[p],"_per_ch_",echo[k],"_common",".nii"))
                    # All active voxels
                    # @infiltrate
                    # tmp_data=niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p],"_girf_",echo[k],"/",scans[p],"_per_ch",".nii"))
                    # Percent signal change
                    # tmp_data=niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p],"_girf_",echo[k],"/psc_msk.nii"))
                    # Delta S over std
                    # tmp_data=niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p],"_girf_",echo[k],"/del_sig_std_msk.nii"))
                    # Delta S over residual std
                    tmp_data=niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p],"_girf_",echo[k],"/del_sig_std_res_bold_msk.nii"))
                    # Delta S
                    # tmp_data=niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p],"_girf_",echo[k],"/del_sig_msk.nii"))

                    tmp_data=tmp_data.raw

                    # Masking only ROI...
                    mask = niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p],"_girf_",echo[k],"/roi_mask.nii"))
                    tmp_data = tmp_data.*mask
                    tmp_data = tmp_data[:]

                    # tmp_data=niread(string("./data/",directories[i],"/analysis/",echos[j],"v_",scans[i],"/tSNR_",metric[k],".nii"))           
                    # tmp_data = tmp_data[tmp_data.!=0]

                    # @infiltrate

                    tmp_data[tmp_data.<0] .= 0
                    tmp_data = tmp_data[tmp_data.!=0]
                    
                    if k==1
                        density!(fig[1,1],tmp_data, color= (:red,0.3), direction =:y, offset = 1.2*p)
                    else
                        density!(fig[1,1],tmp_data, color= (:blue,0.3), direction =:y, offset = 1.2*(p-1))
                    end

                    @info(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p][1:6],"_te_comparison","/",scans[p],"_per_ch_",echo[k],"_common",".nii"))
                    @info(string("Mean: ", mean(tmp_data)))
                    save("./data/tmp/sample.eps",fig)
                    @infiltrate
            end
        # end
    end
end

# @infiltrate
# # Percent signal change
# save(string("./data/tmp/multi_te_perc_actie_vox.eps"),fig)

# Delta s over tSNR
save(string("./data/tmp/multi_te_del_sig_std_res.eps"),fig)

# @infiltrate
# fig = Figure(size = (1600,400))
# ax = Axis(fig[1,1], yticklabelsize=32,
#     limits = (nothing,(0,30)))
#     # ylabel="z-score", ylabelsize=40)    # Only positive
#     # Label(fig[1,1], "sample")
# hidexdecorations!(ax, grid = false)

# violin!(category,data, dodge=dodge, side=side, show_mean=true, color=color)

# save(string("./data/tmp/",scans[i],"_",metric[k],"_raincloud.eps"),fig)
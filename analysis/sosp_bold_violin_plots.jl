using Pkg

cd("/neurodesktop-storage/5T4/Alejandro/sosp_vaso/")
Pkg.activate("./analysis/")

using NIfTI
using CairoMakie
using MIRTjim: jim, prompt; jim(:prompt, true)
using Infiltrator
using Colors
using Statistics

# List all directories.. it will plot 1-2 per each two-sided violin plot
scans = [
        "sb_311_DS_SO_06mm_18fovz_12te_6te_girf_ech2",
        "sb_311_DS_SO_06mm_18fovz_12te_6te_girf_ech1",
        "sb_312_DS_SO_06mm_18fovz_12te_6te_girf_ech2",
        "sb_312_DS_SO_06mm_18fovz_12te_6te_girf_ech1",
        "sb_121_DS_SO_06mm_18fovz_12te_6te_girf_ech2",
        "sb_121_DS_SO_06mm_18fovz_12te_6te_girf_ech1",
        "sb_122_DS_SO_06mm_18fovz_12te_6te_girf_ech2",
        "sb_122_DS_SO_06mm_18fovz_12te_6te_girf_ech1",
        "sb_601_DS_SO_06mm_18fovz_12te_6te_girf_ech2",
        "sb_601_DS_SO_06mm_18fovz_12te_6te_girf_ech1",
        ]

directories = [
                "02062025_sb_9T_paper",
                "02062025_sb_9T_paper",
                "02132025_sb_9T_paper",
                "02132025_sb_9T_paper",
                "02202025_sb_9T_paper",
                ]

drive = ["5T4","5T4","5T4","5T4","5T4"]

categories = 1:0.01:10
# echos=["ech1","ech2"]requ
metric=["tSNR_ROI","clustered_BOLD_ROI"]       # tSNR_ROI, clustered_BOLD_ROI
z_threshold = 1.6

colors=[RGBA{Float32}(1,0,0,0.3),RGBA{Float32}(0,0,1,0.3)]
sides=[:left,:right]

plot_range=[30,10]  # [30,50] for tSNR, []


#####
# tSNR distribution in GM over subjects (Violin plots)
for k=1:length(metric)
p=0
####### BOLD/VASO
data = Vector{Float32}(undef,0)
category =  Vector{Float32}(undef,0)
color = Vector{RGBA{Float32}}(undef,0)
dodge = Vector{Int64}(undef,0)
side = Vector{Symbol}(undef,0)
    for i=1:length(directories)
        for j=1:2  # Two sided Violing plots...

            p += 1

            # Spiral, RED
            tmp_data=niread(string("/neurodesktop-storage/",drive[i],"/Alejandro/sosp_vaso/data/",directories[i],"/analysis/",scans[p],"/",metric[k],".nii"))
            # tmp_data=niread(string("./data/",directories[i],"/analysis/",echos[j],"v_",scans[i],"/tSNR_",metric[k],".nii"))           
            tmp_data=tmp_data.raw
            # msk=niread(string("./data/",directories[i],"/analysis/",echos[j],"v_",scans[i],"/","mask.nii"))
            # msk=msk.raw

            # tmp_data=tmp_data.*msk

            # gm_roi_msk=niread(string("./data/",directories[i],"/analysis/",echos[j],"v_",scans[i],"/",echos[j],"v_",scans[i],"_gm_msk.nii"))
            # gm_roi_msk=gm_roi_msk.raw
            # gm_roi_msk[gm_roi_msk.<50].=0
            # gm_roi_msk[gm_roi_msk.>0].=1
            # tmp_data=tmp_data.*gm_roi_msk

            # Threshold z-scores...
            if metric[k] == "clustered_BOLD_ROI"
                tmp_data[tmp_data.<z_threshold].=0
            end

            tmp_data[tmp_data.>1000].=0
            tmp_data=tmp_data[tmp_data.>0]

            tmp_data=convert(Array{Float32},tmp_data)
            # println("Mean: ",echos[j],"v_",scans[i], " ", metric[k],"   ",mean(tmp_data[:]))
            append!(data,tmp_data)

            tmp_data_category=Vector{Float32}(undef,length(tmp_data))
            tmp_data_category.=categories[i]
            append!(category,tmp_data_category)

            tmp_color =  Vector{RGBA{Float32}}(undef,length(tmp_data))
            tmp_color .= colors[j]
            append!(color,tmp_color)

            tmp_dodge = Vector{Int16}(undef,length(tmp_data))
            tmp_dodge .= i
            append!(dodge,tmp_dodge)

            tmp_side = Vector{Symbol}(undef,length(tmp_data))
            tmp_side .= sides[j]
            append!(side,tmp_side)
            
        end
    end

# fig = Figure(size = (1600,800))         # For 5 "subjects"
fig = Figure(size = (1200,500))         # For 3 "subjects"
ax = Axis(fig[1,1], yticklabelsize=64,
    limits = (nothing,(0,plot_range[k])))
    # ylabel="z-score", ylabelsize=40)    # Only positive
    # Label(fig[1,1], "sample")
hidexdecorations!(ax, grid = false)

violin!(category,data, dodge=dodge, side=side, color=color)

# @infiltrate

save(string("./data/tmp/",scans[1],"_",metric[k],"_violin.eps"),fig)

end

# @infiltrate
# fig = Figure(size = (1600,400))
# ax = Axis(fig[1,1], yticklabelsize=32,
#     limits = (nothing,(0,30)))
#     # ylabel="z-score", ylabelsize=40)    # Only positive
#     # Label(fig[1,1], "sample")
# hidexdecorations!(ax, grid = false)

# violin!(category,data, dodge=dodge, side=side, show_mean=true, color=color)

# save(string("./data/tmp/",scans[i],"_",metric[k],"_raincloud.eps"),fig)
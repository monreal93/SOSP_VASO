using Pkg

cd("/usr/share/5T4/Alejandro/sosp_vaso/")
Pkg.activate("./analysis/")

using NIfTI
using CairoMakie
using MIRTjim: jim, prompt; jim(:prompt, true)
using Infiltrator
using Colors
using Statistics

scans = ["01","41","51","71","81"]
directories = ["05192023_sv_paper","05252023_sv_paper","08232023_sv_abc_paper","09192023_sv_abc_paper","10252023_sv_paper"]
categories = 1:0.01:10
contrast=["v","b"]  # BOLD, VASO
readouts=["s","c"]
colors=[RGBA{Float32}(1,0,0,0.3),RGBA{Float32}(0,0,1,0.3)]
sides=[:left,:right]

plot_range=[30,50]  # 25, 15


#####
# tSNR distribution in GM over subjects (Violin plots)

for k=1:length(contrast)
####### BOLD/VASO
data = Vector{Float32}(undef,0)
category =  Vector{Float32}(undef,0)
color = Vector{RGBA{Float32}}(undef,0)
dodge = Vector{Int64}(undef,0)
side = Vector{Symbol}(undef,0)
    for i=1:length(scans)
        for j=1:length(readouts)


            # Spiral, RED
            tmp_data=niread(string("./data/",directories[i],"/analysis/",readouts[j],"v_",scans[i],"/eff_tSNR_",contrast[k],".nii"))
            # tmp_data=niread(string("./data/",directories[i],"/analysis/",readouts[j],"v_",scans[i],"/tSNR_",contrast[k],".nii"))           
            tmp_data=tmp_data.raw
            msk=niread(string("./data/",directories[i],"/analysis/",readouts[j],"v_",scans[i],"/","mask.nii"))
            msk=msk.raw

            tmp_data=tmp_data.*msk

            # gm_roi_msk=niread(string("./data/",directories[i],"/analysis/",readouts[j],"v_",scans[i],"/",readouts[j],"v_",scans[i],"_gm_msk.nii"))
            # gm_roi_msk=gm_roi_msk.raw
            # gm_roi_msk[gm_roi_msk.<50].=0
            # gm_roi_msk[gm_roi_msk.>0].=1
            # tmp_data=tmp_data.*gm_roi_msk

            tmp_data[tmp_data.>1000].=0
            tmp_data=tmp_data[tmp_data.>0]
            tmp_data=convert(Array{Float32},tmp_data)
            println("Mean: ",readouts[j],"v_",scans[i], " ", contrast[k],"   ",mean(tmp_data[:]))
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

fig = Figure(size = (1600,800))
ax = Axis(fig[1,1], yticklabelsize=64,
    limits = (nothing,(0,plot_range[k])))
    # ylabel="z-score", ylabelsize=40)    # Only positive
    # Label(fig[1,1], "sample")
hidexdecorations!(ax, grid = false)

violin!(category,data, dodge=dodge, side=side, show_mean=true, color=color)

@infiltrate

save(string("./data/tmp/",contrast[k],"_tSNR_violin.eps"),fig)

end

# @infiltrate
# fig = Figure(size = (1600,400))
# ax = Axis(fig[1,1], yticklabelsize=32,
#     limits = (nothing,(0,30)))
#     # ylabel="z-score", ylabelsize=40)    # Only positive
#     # Label(fig[1,1], "sample")
# hidexdecorations!(ax, grid = false)

# violin!(category,data, dodge=dodge, side=side, show_mean=true, color=color)

# save(string("./data/tmp/",scans[i],"_",contrast[k],"_raincloud.eps"),fig)
using Pkg

cd("/usr/share/5T4/Alejandro/sosp_vaso/")
Pkg.activate("./analysis/")

using NIfTI
using CairoMakie
using MIRTjim: jim, prompt; jim(:prompt, true)
using Infiltrator
using Colors

scans = ["01","41","51","71","81"]
directories = ["05192023_sv_paper","05252023_sv_paper","08232023_sv_abc_paper","09192023_sv_abc_paper","10252023_sv_paper"]
categories = ["A","B","C","D","E","F","G","H","I","J","K"]
contrast="VASO"  # BOLD, VASO

if contrast=="VASO"
    plot_range=15  # 25, 15
    plot_threshold = 2.3
else
    plot_range=25  # 25, 15
    plot_threshold = 2.3
end
####### BOLD/VASO
data = Vector{Float32}(undef,0)
category = Vector{String}(undef,0)
color = Vector{RGBA{Float32}}(undef,0)
dodge = Vector{Int64}(undef,0)
j=1
for i=1:length(scans)

    fig = Figure(size = (400,800))
    ax = Axis(fig[1,1], yticklabelsize=32,
        ygridvisible = false,
        # limits = (nothing,(-plot_range,plot_range)))  # Including negative values
        limits = (nothing,(0,plot_range)))
        # ylabel="z-score", ylabelsize=40)    # Only positive
    # Label(fig[1,1], "sample")
    hidexdecorations!(ax, grid = false)

    # Spiral, RED
    # tmp_data=niread(string("./data/",directories[i],"/analysis/sv_",scans[i],"/",contrast,"_msk.nii")) # raw maps
    tmp_data=niread(string("./data/",directories[i],"/analysis/sv_",scans[i],"/clustered_",contrast,".nii")) # clustered maps
    per_spi = length(tmp_data)
    # tmp_data=tmp_data[tmp_data.!=0]  # All values different than 0
    tmp_data=tmp_data[tmp_data.>=plot_threshold]  # Thresholded values
    # tmp_data=tmp_data[tmp_data.>0]  # Thresholded values >0
    per_spi = length(tmp_data)/per_spi*100
    # append!(data,tmp_data)
    println(string("Total Spiral pixels ",contrast," ", scans[i],":\t", length(tmp_data),"\t % ",per_spi))

    tmp_data_category=Vector{String}(undef,length(tmp_data))
    tmp_data_category.=categories[i]
    # append!(category,tmp_data_category)
    # j=j+1

    tmp_color =  Vector{RGBA{Float32}}(undef,length(tmp_data))
    tmp_color .= RGBA{Float32}(1,0,0,0.3)
    # append!(color,tmp_color)

    tmp_dodge = Vector{Int16}(undef,length(tmp_data))
    tmp_dodge .= i
    # append!(dodge,tmp_dodge)

    rainclouds!(tmp_data_category,tmp_data[:],cloud_width=1, dodge=tmp_dodge, color=tmp_color)
    
    # EPI, BLUE
    # tmp_data=niread(string("./data/",directories[i],"/analysis/cv_",scans[i],"/",contrast,"_msk.nii"))  # raw maps
    tmp_data=niread(string("./data/",directories[i],"/analysis/cv_",scans[i],"/clustered_",contrast,".nii")) # clustered maps
    per_epi = length(tmp_data)
    # tmp_data=tmp_data[tmp_data.!=0]  # All values different than 0
    tmp_data=tmp_data[tmp_data.>=plot_threshold]  # Thresholded values
    # tmp_data=tmp_data[tmp_data.>0]  # Thresholded values >0
    per_epi = length(tmp_data)/per_epi*100
    # append!(data,tmp_data)
    println(string("Total Cartesian pixels ",contrast," ", scans[i],":\t", length(tmp_data),"\t % ",per_epi))

    tmp_data_category=Vector{String}(undef,length(tmp_data))
    tmp_data_category.=categories[i]
    # append!(category,tmp_data_category)
    # j=j+1

    tmp_color =  Vector{RGBA{Float32}}(undef,length(tmp_data))
    tmp_color .= RGBA{Float32}(0,0,1,0.3)
    # append!(color,tmp_color)

    tmp_dodge = Vector{Int16}(undef,length(tmp_data))
    tmp_dodge .= i
    # append!(dodge,tmp_dodge)
    
    # Plotting the EPI pdf compared to spiral pdf (% of active voxels)
    cloud_width = per_epi/per_spi
    rainclouds!(tmp_data_category,tmp_data[:], cloud_width=cloud_width, dodge=tmp_dodge, color=tmp_color)

    @infiltrate

    save(string("./data/tmp/",scans[i],"_",contrast,"_raincloud.eps"),fig)
end
#########


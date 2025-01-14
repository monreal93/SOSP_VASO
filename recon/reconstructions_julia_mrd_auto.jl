using Pkg

cd("/usr/share/5T4/Alejandro/sosp_vaso/")
# cd("/usr/share/5T3/Alejandro/sosp_vaso/")
Pkg.activate("./recon/")

using Revise
# using Infiltrator
using MRIReco, MRIFiles, MRICoilSensitivities, MRIFieldmaps
using MAT, NIfTI , JLD #, MriResearchTools
using StatsBase: mean
using Unitful: s, ustrip
using MIRTjim: jim, prompt; jim(:prompt, true)
using ImageTransformations, Interpolations
using NFFT
using FLoops
using Statistics
using DSP
# using ROMEO

using FFTW
FFTW.set_provider!("mkl")
# FFTW.set_provider!("fftw")

include("./functions/fn_save_nii.jl")
include("./functions/fn_Utilities.jl")
include("./functions/fn_PrepareMRD.jl")
include("./functions/fn_ReconstructCartesian.jl")
include("./functions/fn_ReconstructSpiral.jl")
include("./functions/fn_RawDataCorrection.jl")
include("./functions/fn_CalculateSensitivityOffresonanceMaps.jl")
# include("./functions/fn_calculateSphericalHarmonics.jl")
# include("./functions/fn_motionCorrection.jl")
include("./functions/fn_CorrectOffResonance.jl")

params = Dict{Symbol, Any}()
params[:do_pi_recon] = true             # Perform PI reconstruction or direct recon
params[:do_b0_corr] = false
params[:do_b0_corr_seg] = 0   # Perform custom  0=Gridding(MRIReco), 1=Time Segmented, 2=Frequency segmented,
params[:do_rDORK_corr] = false      
params[:do_iDORK_corr] = false            # Perform Interleaf DORK (WIP)
params[:do_coco_corr] = false           # Perform concomitant field correction (WIP)  
params[:do_k0_corr] = false             # Perform K0 demodulation of raw data, only works with sk/girf trajectory (WIP)
params[:do_high_b0_corr] = false         # Higher order correction, FID navigators (WIP) 
params[:do_motion_corr] = false         # Motion correction, FID navigators (WIP)
params[:do_t2s_corr] = false
params[:is2d] = false
params[:rep_recon] = 4            # Range of rep to recon, 0 for all rep, ex. (5:5)- rep 5 
params[:traj_type] = "nom"                 # Trajectory type nominal="nom",skope="sk",poet = "poet", girf = "girf"
params[:save_ph] = false                       # Save phase of recon as nifti
params[:mcorr] = ""           # Motion correction with navigators "_mCorr"
params[:recon_order] =  1                  # Higher order recon (2,3)

# Some parameters
# params[:scan] = "sb_004_DS_SO_08mm_2te"            # For now: if multipe echos, include _e1.. _e2..
params[:scan] = ARGS[2]
params[:scan_b0] = "s01"           # Name of the ME-GRE to use for CS and B0map, ONLY USED WHEN SEPARATE FIELDMAP SCAN...
# params[:directory] = "10282024_sb_9T"        # directory where the data is stored'
params[:directory] = ARGS[1]

# Drive where data is stored... (5T4/5T3)
path_tmp = "/usr/share/5T3/Alejandro/sosp_vaso/"
params[:path] = string(path_tmp,"data/",params[:directory])

# Read pulseq parameters
params_pulseq = matread(string(params[:path],"/acq/",params[:scan],"_params.mat")); params_pulseq = params_pulseq["params"]

# Some constraints
if params[:recon_order] > 1; params[:traj_type] = "sk"; end
# rDORK only for spiral OUT
# if params_pulseq["spi"]["type"] != 0; params[:do_rDORK_corr] = false; end

# File names
if params_pulseq["gen"]["me_gre"] > 0
    idx = findall.("_",params[:scan])[1][1]+1:findall.("_",params[:scan])[2][1]-1
    file_name = string(string(params[:scan][1],params[:scan][idx]),"_",Int(params_pulseq["gen"]["n_ov"][1]),"_",Int(params_pulseq["gen"]["n_ov"][2]),"_",Int(params_pulseq["gen"]["n_ov"][3]))
else
    file_name = string(params[:scan_b0],"_",Int(params_pulseq["gen"]["n_ov"][1]),"_",Int(params_pulseq["gen"]["n_ov"][2]),"_",Int(params_pulseq["gen"]["n_ov"][3]))
end
save_suffix = string("_",params[:traj_type]) 

# read spiral raw MRD file, add partition labels, format into expected array shape
file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/",params[:scan],".h5"))
@info("Reading raw data...")
rawData = RawAcquisitionData(file)

# Splitting raw data into fieldmap and fMRI, and calculate Drift
RepRef = 1         # DORK/OffresonanceDrift reference repetition
B0RepRef = 1       # If multiple Repetitions acquired of ME-GRE, choose wihich one to use for Drift and Sensitivity/Fieldmap
if params_pulseq["gen"]["me_gre"] == 1
    # Read pulseq parameters
    params_me_gre_pulseq = matread(string(params[:path],"/acq/",params[:scan],"_params_me_gre.mat")); params_me_gre_pulseq = params_me_gre_pulseq["params_me_gre"]
    MeGreProfiles = Int(params_me_gre_pulseq["gen"]["n_ov"][3]*params_me_gre_pulseq["spi"]["interl"]*params_me_gre_pulseq["gen"]["me_gre_echos"])

    rawData_b0 = RawAcquisitionData(rawData.params,rawData.profiles[1:MeGreProfiles])
    rawData_b0.params["TE"] = vec(params_me_gre_pulseq["gen"]["te"].*1e3)

    # Remove the profiles from ME-GRE
    rawData.profiles = rawData.profiles[Int(MeGreProfiles+1):end]
    
    numInterl_b0 = Int(params_me_gre_pulseq["spi"]["interl"])
    numInterl = Int(params_pulseq["spi"]["interl"])
    numPar = Int(params_me_gre_pulseq["gen"]["n_ov"][3]./params_me_gre_pulseq["gen"]["kz"])
    numEch = Int(params_me_gre_pulseq["gen"]["me_gre_echos"])
    # Sets of functional data
    numSet = maximum(unique([rawData.profiles[l].head.idx.set+1 for l=1:length(rawData.profiles)]))

    if params_pulseq["spi"]["type"] == 0
        # Drift from fieldmap to first acquisition, First mean then angle...
        # nav = mean(rawData.profiles[Int(((numPar/2*nkumSet)+1))].data[6:14,:], dims=2)  # Do I want to get te second or first repetition?
        nav = mean(rawData.profiles[Int(((RepRef-1)*numPar*numSet)+(((numPar/2+1)*numSet*numInterl)+1))].data[3:8,:], dims=2)  # Do I want to get te second or first repetition?
        nav_b0 = mean(rawData_b0.profiles[Int((((numPar/2+1)*numInterl_b0*numEch)+1))].data[3:8,:], dims=2)
        nav_m = mean(unwrap(angle.(nav); dims=1))
        nav_b0_m = mean(unwrap(angle.(nav_b0); dims=1))

        OffResonanceDrift = (nav_m-nav_b0_m)./params_me_gre_pulseq["gen"]["te"][1]
    end

elseif params_pulseq["gen"]["me_gre"] == 2
    file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/",params[:scan],"_me_gre.h5"))
    rawData_b0 = RawAcquisitionData(file)

    # Read pulseq parameters
    params_me_gre_pulseq = matread(string(params[:path],"/acq/",params[:scan],"_params_me_gre.mat")); params_me_gre_pulseq = params_me_gre_pulseq["params_me_gre"]
    MeGreProfiles = Int(params_me_gre_pulseq["gen"]["n_ov"][3]*params_me_gre_pulseq["spi"]["interl"]*params_me_gre_pulseq["gen"]["me_gre_echos"])

    numInterl_b0 = Int(params_me_gre_pulseq["spi"]["interl"])
    numInterl = Int(params_pulseq["spi"]["interl"])
    numPar = Int(params_me_gre_pulseq["gen"]["n_ov"][3]./params_me_gre_pulseq["gen"]["kz"])
    numEch = Int(params_me_gre_pulseq["gen"]["me_gre_echos"])
    # Sets of functional data
    numSet = maximum(unique([rawData.profiles[l].head.idx.set+1 for l=1:length(rawData.profiles)]))

    if params_pulseq["spi"]["type"] == 0
        # Drift from fieldmap to first acquisition, First mean then angle...
        # nav = mean(rawData.profiles[Int(((numPar/2*numSet)+1))].data[6:14,:], dims=2)  # Do I want to get te second or first repetition?
        nav = mean(rawData.profiles[Int(((RepRef-1)*numPar*numSet)+(((numPar/2+1)*numSet*numInterl)+1))].data[3:8,:], dims=2)  # Do I want to get te second or first repetition?
        # If multiple rep of ME-GRE, use B0RepRef to select the desired one
        nav_b0_profile = Int((((numPar/2+1)*numInterl_b0*numEch)+1))
        nav_b0_profile = nav_b0_profile .+ (MeGreProfiles*(B0RepRef-1))
        nav_b0 = mean(rawData_b0.profiles[nav_b0_profile].data[3:8,:], dims=2)
        nav_m = mean(unwrap(angle.(nav); dims=1))
        nav_b0_m = mean(unwrap(angle.(nav_b0); dims=1))

        OffResonanceDrift = (nav_m-nav_b0_m)./params_me_gre_pulseq["gen"]["te"][1]
    end
end

if isfile(string(params[:path],"/acq/me_gre_",file_name,".jld"))
    @info ("ME GRE exists... Reconstruct again? y/n")
    input = readline()
else
    input = ""
end
if input == "n"
    @info ("Loading ME GRE ...")
    recon_b0 = load(string(params[:path],"/acq/me_gre_",file_name,".jld"))
    recon_b0 = recon_b0["recon_b0"]
    if params_pulseq["gen"]["me_gre"] != 0
        params[:b0_TE] = vec(params_me_gre_pulseq["gen"]["te"].*1e3)
    end
else
    @info ("Reconstructing ME GRE ...")
    # Reconstruction of ME-GRE scan (spiral from pulseq or cartesian separate scan)
    if params_pulseq["gen"]["me_gre"] > 0
        ks_traj_me_gre = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_me_gre_",params[:traj_type],".mat")); ks_traj_me_gre = ks_traj_me_gre["ks_traj_me_gre"]
        # Normalizing
        # Not sure if I want to do -1 for kx
        ks_traj_me_gre["kx"] = ks_traj_me_gre["kx"]./(maximum([abs(minimum(ks_traj_me_gre["kx"])),abs(maximum(ks_traj_me_gre["kx"]))])*2)
        ks_traj_me_gre["ky"] = ks_traj_me_gre["ky"]./(maximum([abs(minimum(ks_traj_me_gre["ky"])),abs(maximum(ks_traj_me_gre["ky"]))])*2) .*(-1)
        ks_traj_me_gre["kz"] = ks_traj_me_gre["kz"]./(maximum([abs(minimum(ks_traj_me_gre["kz"])),abs(maximum(ks_traj_me_gre["kz"]))])*2)
        ks_traj_me_gre = hcat(ks_traj_me_gre["kx"][:],ks_traj_me_gre["ky"][:],ks_traj_me_gre["kz"][:])
        ks_traj_me_gre = permutedims(ks_traj_me_gre,[2,1])

        ks_traj_me_gre = convert(AbstractMatrix{Float32},ks_traj_me_gre)
        times_me_gre = repeat(params_me_gre_pulseq["gen"]["t_vector"][:],Int(params_me_gre_pulseq["gen"]["n_ov"][3]))
        times_me_gre = convert(Vector{Float32},times_me_gre)

        ks_traj_me_gre = Trajectory(ks_traj_me_gre,1,Int(params_me_gre_pulseq["gen"]["ro_samples"]*params_me_gre_pulseq["spi"]["interl"]); 
                times=times_me_gre,TE=params_me_gre_pulseq["gen"]["TE"],AQ=params_me_gre_pulseq["gen"]["ro_time"], 
                numSlices=Int(params_me_gre_pulseq["gen"]["n_ov"][3]),circular=true)
        
        # Reconstruct echo by echo
        kdata = Array{Array{ComplexF32,2},3}(undef,1,1,1)
        recon_b0 = Array{ComplexF32,5}(undef,(Int.(params_me_gre_pulseq["gen"]["n_ov"])...,size(rawData.profiles[1].data,2),Int.(params_me_gre_pulseq["gen"]["me_gre_echos"])))
        for i_ech=1:params_me_gre_pulseq["gen"]["me_gre_echos"]

            b0_ech_range = Int(i_ech):Int(params_me_gre_pulseq["gen"]["me_gre_echos"]):MeGreProfiles
            b0_ech_range = b0_ech_range .+ (MeGreProfiles*(B0RepRef-1))

            rawData_b0_ech = RawAcquisitionData(rawData_b0.params,rawData_b0.profiles[b0_ech_range])

            # rawData_b0 = RawAcquisitionData(rawData.params,rawData.profiles[Int(i_ech):Int(params_me_gre_pulseq["gen"]["me_gre_echos"]):MeGreProfiles])
            rawData_b0_ech.params["TE"] = vec(params_me_gre_pulseq["gen"]["te"].*1e3)

            numRead = size(rawData_b0_ech.profiles[1].data,1) # Readout from raw data
            numInterl = Int(params_me_gre_pulseq["spi"]["interl"])
            numPar = Int(params_me_gre_pulseq["gen"]["n_ov"][3]./params_me_gre_pulseq["gen"]["kz"])
            numSet = maximum(unique([rawData_b0_ech.profiles[l].head.idx.set+1 for l=1:length(rawData_b0_ech.profiles)]))
            numCha = size(rawData_b0_ech.profiles[1].data,2)
            numRep = rawData_b0_ech.params["userParameters"]["sWipMemBlock.alFree[7]"]  # I think this is repetitions in special card
            # Adding some extra parameters
            params[:numRead] = numRead*numSet
            params[:numInterl] = numInterl
            params[:numSet] = numSet
            params[:numPar] = numPar
            params[:numCha] = numCha
            params[:numRep] = numRep

            tmp = FormatRawData(rawData_b0_ech,params;single_rep=true,rep_format=1)

            kdata[1] = dropdims(reshape(tmp,(:,numCha,1)),dims=3)
            
            acqData_b0 = AcquisitionData(ks_traj_me_gre,kdata; encodingSize=(Int(params_me_gre_pulseq["gen"]["ro_samples"]*params_me_gre_pulseq["spi"]["interl"]),1,Int(params_me_gre_pulseq["gen"]["n_ov"][3])), fov=Tuple(params_me_gre_pulseq["gen"]["fov"]))

            recon_b0[:,:,:,:,Int(i_ech)] = ReconSpiralMEGRE(acqData_b0,params_me_gre_pulseq)

        end
        params[:b0_TE] = vec(params_me_gre_pulseq["gen"]["te"].*1e3)
    else
        # read ME-GRE MRD data and recon
        file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/b0_",params[:scan_b0],"_fieldmap.h5"))
        acqData_b0 = AcquisitionData(file)
        rawData_b0 = RawAcquisitionData(file)

        if params_pulseq["gen"]["field_strength"] == 7im
            recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=false)
        elseif params_pulseq["gen"]["field_strength"] == 7
            recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=true)
        elseif params_pulseq["gen"]["field_strength"] == 9
            recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=true) 
        end
    end
    save(string(params[:path],"/acq/me_gre_",file_name,".jld"),"recon_b0",recon_b0)

    if params_pulseq["gen"]["me_gre"] == 0
        # # AMM: For now I add them manually... Fix it..
        @info("Manually write the TE for B0...")
        params[:b0_TE] = vec([2.99,6.31,12,25,32])
    end

    # Saving ME GRE for reference in tmp folder
    b0_nii = recon_b0[:,:,:,:,:]
    b0_nii = dropdims(sqrt.(sum(abs.(b0_nii).^2, dims=4)), dims=4)
    b0_nii = imresize(b0_nii,Int(params_pulseq["gen"]["n"][1]),Int(params_pulseq["gen"]["n"][2]),Int(params_pulseq["gen"]["n"][3]))
    niwrite(string(params[:path],"/tmp/",params[:scan],".nii"),NIVolume(b0_nii; voxel_size=Tuple(params_pulseq["gen"]["res"].*1e3)))
end
recon_b0_sos = dropdims(sqrt.(((sum(abs.(recon_b0).^2; dims=4)))), dims=4)

# Getting some parameters from the raw data
numRead = size(rawData.profiles[1].data,1) # Readout from raw data
# numRead = Int(params_pulseq["gen"]["ro_samples"])  # Readout from pulseq params
numInterl = Int(params_pulseq["spi"]["interl"])
# numPar = maximum(unique([rawData.profiles[l].head.idx.kspace_encode_step_2+1 for l=1:length(rawData.profiles)]))
numPar = Int(params_pulseq["gen"]["n_ov"][3]./params_pulseq["gen"]["kz"])
numSet = maximum(unique([rawData.profiles[l].head.idx.set+1 for l=1:length(rawData.profiles)]))
numCha = size(rawData.profiles[1].data,2)
numRep = rawData.params["userParameters"]["sWipMemBlock.alFree[7]"]  # I think this is repetitions in special card
# numRep = Int(length(rawData.profiles)/numPar/numSet)

# Adding some extra parameters
params[:numRead] = numRead*numSet
params[:numInterl] = numInterl
params[:numSet] = numSet
params[:numPar] = numPar
params[:numCha] = numCha
params[:numRep] = numRep
params[:kz_enc] = params_pulseq["gen"]["kz_enc"]
params[:TE] = params_pulseq["gen"]["TE"]
params[:acq_times] = params_pulseq["gen"]["t_vector"].+6e-3
params[:reconSize] = (params_pulseq["gen"]["n_ov"][1],params_pulseq["gen"]["n_ov"][2],params_pulseq["gen"]["n_ov"][3])
params[:FieldOfView] = params_pulseq["gen"]["fov"]


if params_pulseq["gen"]["seq"] == 1
    # SS-SI-VASO
    params[:contrasts] = ["v","b"]
    params[:numRep] = params[:numRep]*2
elseif params_pulseq["gen"]["seq"] == 2
    # ABC
    params[:contrasts] = ["abc"]
elseif params_pulseq["gen"]["seq"] == 4
    # BOLD
    params[:contrasts] = ["b"]
end

# Calculate Sensitivity maps
if params[:do_pi_recon]
    save_suffix = string(save_suffix,"_cs")
    if isfile(string(params[:path],"/acq/cs_",file_name,".jld"))
        @info ("Sensitivity map exists... Re-calculate? y/n")
        input = readline()
    else
        input = ""
    end
    if input == "n"
        @info ("Loading Sensitivity maps ...")
        SensitivityMap = load(string(params[:path],"/acq/cs_",file_name,".jld"))
        SensitivityMap = SensitivityMap["SensitivityMap"]
    else
        @info("Calculating Sensitivity Maps...")
        SensitivityMap = CalculateSensitivityMap(recon_b0,Tuple(Int.(params_pulseq["gen"]["n_ov"])))
        save(string(params[:path],"/acq/cs_",file_name,".jld"),"SensitivityMap",SensitivityMap)
    end
end

# Calculate B0 map
if params[:do_b0_corr] || (params[:traj_type] == "sk" && params[:recon_order] > 1)
    if params[:do_high_b0_corr] && params[:do_b0_corr_seg] == 0
        save_suffix = string(save_suffix,"_hb0")
    elseif params[:do_b0_corr_seg] == 1
        save_suffix = string(save_suffix,"_tsb0")
    elseif params[:do_b0_corr_seg] == 2
        save_suffix = string(save_suffix,"_fsb0")
    else
        save_suffix = string(save_suffix,"_b0")
    end
    if isfile(string(params[:path],"/acq/b0_",file_name,".jld"))
        @info ("OffResonance map exists... Re-calculate? y/n")
        input = readline()
    else
        input = ""
    end
    if input == "n"
        @info("Loading Offresonance Maps....")
        OffResonanceMap = load(string(params[:path],"/acq/b0_",file_name,".jld"))
        OffResonanceMap = OffResonanceMap["OffResonanceMap"]
    else
        @info("Calculating Offresonance Maps....")
        OffResonanceMap = CalculateOffresonanceMap(recon_b0,SensitivityMap,params[:b0_TE])
        save(string(params[:path],"/acq/b0_",file_name,".jld"),"OffResonanceMap",OffResonanceMap)
    end
end

# Calculate concomitant field
if params[:do_coco_corr] && params[:do_b0_corr]
    @info("Calculating Concomitant Maps....")
    save_suffix = string(save_suffix,"_co")
    # phase,read,slice
    RotMatrix =vcat(transpose(collect(rawData_b0.profiles[1].head.read_dir)),
            transpose(collect(rawData_b0.profiles[1].head.phase_dir)),
            transpose(collect(rawData_b0.profiles[1].head.slice_dir)))
    CenterPosition = (rawData_b0.profiles[end].head.position .+ rawData_b0.profiles[1].head.position)./2 .* 1e-3
    CocoFieldMap = CalculateConcomitantFieldMap(RotMatrix,CenterPosition, params_pulseq)
    # CocoFieldMap = CocoFieldMap .* (2*π)
    # CocoFieldMap = CocoFieldMap .- Float32(10*2π) # Ading some extra Hz, maybe due to drift?
    # I need to reverse dims=3,2,1, maybe rot matrix is wrong?
    # CocoFieldMap = reverse(CocoFieldMap,dims=(1,2,3)) # maybe just for 7T?
    # OffResonanceMap .= OffResonanceMap .* (π/2)
    # Sometimes - works better, sometimes +
    OffResonanceMap .= OffResonanceMap .- (CocoFieldMap.*im)
    # OffResonanceMap .= OffResonanceMap .- (14*2*π*im)
    mask = SensitivityMap[:,:,:,1]
    mask[abs.(mask) .> 0] .= 1
    mask = isone.(mask)
    OffResonanceMap .= OffResonanceMap .* mask
end

# Higher order B0 correction, FID navigators approach
if params[:do_high_b0_corr] && params[:do_b0_corr]
    @info("Calculating δB0 Maps....")
    # Calibration matrix
    A, B, sh_basis, ΔB0, b_A = getCalibrationMatrix(params,OffResonanceMap,recon_b0)
    # @info("Stop... high B0 corr...")
    # @infiltrate
end

# Normalize k-space trajectory and create Trajectory object
if params[:traj_type] == "sk"
    ks_traj = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_",params[:traj_type],".mat")); ks_traj = ks_traj["ks_traj"]
    # ToDo: Sync Skope trajectory to raw data...

    ## Approach 2) same structure as nominal
    # Normalizing
    ks_traj["kx"] = ks_traj["kx"]./(maximum([abs(minimum(ks_traj["kx"])),abs(maximum(ks_traj["kx"]))])*2)
    ks_traj["ky"] = ks_traj["ky"]./(maximum([abs(minimum(ks_traj["ky"])),abs(maximum(ks_traj["ky"]))])*2).*(-1)
    ks_traj["kz"] = ks_traj["kz"]./(maximum([abs(minimum(ks_traj["kz"])),abs(maximum(ks_traj["kz"]))])*2)
    # k0 = ks_traj["k0"][:]
    ks_traj = hcat(ks_traj["kx"][:],ks_traj["ky"][:],ks_traj["kz"][:])
    ks_traj = permutedims(ks_traj,[2,1])

else
    ks_traj = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_",params[:traj_type],".mat")); ks_traj = ks_traj["ks_traj"]
    # Normalizing
    ks_traj["kx"] = ks_traj["kx"]./(maximum([abs(minimum(ks_traj["kx"])),abs(maximum(ks_traj["kx"]))])*2)
    ks_traj["ky"] = ks_traj["ky"]./(maximum([abs(minimum(ks_traj["ky"])),abs(maximum(ks_traj["ky"]))])*2).*(-1)
    ks_traj["kz"] = ks_traj["kz"]./(maximum([abs(minimum(ks_traj["kz"])),abs(maximum(ks_traj["kz"]))])*2)
    ks_traj = hcat(ks_traj["kx"][:],ks_traj["ky"][:],ks_traj["kz"][:])
    ks_traj = permutedims(ks_traj,[2,1])
end

# K0 from skope
if params[:do_k0_corr] && params[:traj_type] == "sk"
    k0_meas = matread(string(params[:path],"/acq/",params[:scan],"_k0_sk.mat")); k0_meas = k0_meas["k0"][:]
    k0_sim = matread(string(params[:path],"/acq/",params[:scan],"_k0_sim.mat")); k0_sim = k0_sim["k0"]
elseif params[:do_k0_corr] && params[:traj_type] == "girf"
    # @info("Stop... K0..")
    # @infiltrate
    k0_meas = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_",params[:traj_type],".mat")); k0_meas = k0_meas["ks_traj"]
    k0_meas = Float64.(k0_meas["k0"])
    k0_meas = k0_meas.*(-1) # Do I need to do -1, I do it for ks_traj..?
    k0_sim = matread(string(params[:path],"/acq/",params[:scan],"_k0_sim.mat")); k0_sim = k0_sim["k0"]
end

ks_traj = convert(AbstractMatrix{Float32},ks_traj)
times = repeat(params_pulseq["gen"]["t_vector"][:],numPar)
times = convert(Vector{Float32},times)

ks_traj = Trajectory(ks_traj,1,Int(params_pulseq["gen"]["ro_samples"]*params_pulseq["spi"]["interl"]); 
        times=times,TE=params_pulseq["gen"]["TE"],AQ=params_pulseq["gen"]["ro_time"], 
        numSlices=Int(params_pulseq["gen"]["n_ov"][3]),circular=true)


# Setting repetition label
rawData = SetRepetitionLabel(rawData,params)

# K0 correction
if params[:do_k0_corr]
    save_suffix = string(save_suffix,"_k0")
end

# Reconstruction parameters
params_recon = Dict{Symbol, Any}()
params_recon = merge(defaultRecoParams(), params_recon)
if params[:do_pi_recon]
    params_recon[:reco] = "multiCoil"
    params_recon[:senseMaps] = SensitivityMap
end
if params[:do_b0_corr]
    # Temp: Scaling fieldmap for testing...
    # OffResonanceMap = OffResonanceMap .* Float32(π)
    # Temp: Adding some Hz (rad) for drift
    # For 9T I need to do *-1 (not sure why..)
    # @info("Stop... OffResonance...Drift..")
    # @infiltrate
    if params_pulseq["spi"]["type"] == 0
        if params_pulseq["gen"]["me_gre"] != 0
            if params_pulseq["gen"]["field_strength"] == 9
                OffResonanceMap = OffResonanceMap .- ComplexF32(abs(OffResonanceDrift).*im)
                # OffResonanceMap = OffResonanceMap .- ComplexF32(300*im) # Maybe I am missing some Hz/rad..
                # OffResonanceMap = OffResonanceMap .* ComplexF32(2) # AMM: this helps with B0 corr, check...
            else
                OffResonanceMap = OffResonanceMap .- ComplexF32(abs(OffResonanceDrift).*im)
            end
        end
    end

    # AMM: Temp: Adding some rad to fieldmap.. test...
    # OffResonanceMap = OffResonanceMap .- ComplexF32(80*im) 
    if params[:do_b0_corr_seg] == 0
        params_recon[:correctionMap] = deepcopy(OffResonanceMap)
    end
end
if params[:recon_order] > 1 || (params[:traj_type] == "sk" && params[:recon_order] > 1)
    params_recon[:reco] = "expandedSignalModel"
    params_recon[:higherOrderTrajectory] = convert(Matrix{Float32},ks_traj_high)
    params_recon[:senseMaps] = SensitivityMap
    params_recon[:correctionMap] = deepcopy(OffResonanceMap)
end

params_pulseq["gen"]["n_ov"] = Int.(params_pulseq["gen"]["n_ov"])
params_recon[:reconSize] = (params_pulseq["gen"]["n_ov"][1],params_pulseq["gen"]["n_ov"][2],params_pulseq["gen"]["n_ov"][3])
params_recon[:regularization] = "L1"
params_recon[:λ] = 1.e-2  # 1.e-2
params_recon[:iterations] = 40 # (10)
params_recon[:solver] = "admm" # cgnr (L2), fista (L1), admm(L1)
params_recon[:method] = "nfft" # nfft, leastsquare

# FFT Shift
fftShift = (rawData.profiles[1].head.position .- rawData_b0.profiles[1].head.position)./Float32(2)

# repetition DORK correction
if params[:do_rDORK_corr]
    save_suffix = string(save_suffix,"_rDORK")
    # RepRef = 3         # DORK reference repetition
    nav_ref = FormatRawData(rawData,params;single_rep=true,rep_format=RepRef)
    if params_pulseq["spi"]["type"] == 0
        global nav_range = 6:20
    elseif params_pulseq["spi"]["type"] == 1
        global nav_range = Int.(params_pulseq["gen"]["ro_samples"]-19:params_pulseq["gen"]["ro_samples"]-5)
    end
    if params[:kz_enc] == 0  # Linear
        nav_ref = nav_ref[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
        nav_ref = mean(nav_ref, dims=2)
    elseif params[:kz_enc] == 1 # Center-out
        nav_ref = nav_ref[nav_range,1,:,1]
        nav_ref = mean(nav_ref, dims=2)
    end
    if params_pulseq["gen"]["seq"] == 1
        RepRef = RepRef+1         # DORK reference repetition
        nav_ref1 = FormatRawData(rawData,params;single_rep=true,rep_format=RepRef)
        if params[:kz_enc] == 0  # Linear
            nav_ref1 = nav_ref1[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
            nav_ref1 = mean(nav_ref1, dims=2)
        elseif params[:kz_enc] == 1 # Center-out
            nav_ref1 = nav_ref1[nav_range,1,:,1]
            nav_ref1 = mean(nav_ref1, dims=2)
        end
    end
end

# interleave DORK correction
if params[:do_iDORK_corr]
    save_suffix = string(save_suffix,"_iDORK")
    if params_pulseq["spi"]["type"] == 0
        global nav_range = 6:20
    elseif params_pulseq["spi"]["type"] == 1
        global nav_range = Int.(params_pulseq["gen"]["ro_samples"]-19:params_pulseq["gen"]["ro_samples"]-5)
    end
end
###### Temp to test interpolation.. in kspace::
save_suffix_orig = deepcopy(save_suffix)
ks_traj_orig = deepcopy(ks_traj)
######

# Reconstruct by repetition....
kdata = Array{Array{ComplexF32,2},3}(undef,1,1,1)
if params[:rep_recon] == 0:0 || params[:rep_recon] == 0 
    params[:rep_recon] = 1:params[:numRep]
end
@time begin
    for i_rep=params[:rep_recon]

        tmp = FormatRawData(rawData,params;single_rep=true,rep_format=i_rep)

        # FFT Shift
        # ToDo: Use correct value for shift... how do I get it?
        tmp = CorrectFFTShift(tmp,fftShift,ks_traj.nodes,params)

        # K0 correction
        if params[:do_k0_corr] 
            tmp = Correctk0(tmp,k0_meas,k0_sim,params)
        end

        # repetition DORK correction
        if params[:do_rDORK_corr]
            if params_pulseq["gen"]["seq"] == 1
                if isodd(i_rep)
                    tmp = CorrectRepetitionDORK(tmp,nav_ref,nav_range,params)
                else
                    tmp = CorrectRepetitionDORK(tmp,nav_ref1,nav_range,params)
                end
            else
                tmp = CorrectRepetitionDORK(tmp,nav_ref,nav_range,params)
            end
        end

        # interleaf DORK correction
        if params[:do_iDORK_corr]
            for i_int=1:numInterl-1
                tmp = CorrectInterleafDORK(tmp,nav_range,i_int,params)
            end
        end

        # Dynamic Off-resonance correction.. FID navigators
        if params[:do_high_b0_corr] || params[:do_motion_corr]

            nav = tmp[6:14,Int(params_pulseq["gen"]["n_ov"][3]/2+1),:,1]
            nav = mean(nav, dims=1)'
            YY = vcat(real.(nav),imag.(nav))
            AA = vcat(real.(A),imag.(A))

            b = AA \ YY
            δB0 = B * (b)
            δB0 = reshape(δB0,Int.(params[:reconSize]))
            # δB0 = reverse(δB0,dims=(3))
            mask = SensitivityMap[:,:,:,1]
            mask[abs.(mask) .> 0] .= 1
            mask = isone.(mask)
            params_recon[:correctionMap] .= ((deepcopy(OffResonanceMap) .+ real.(δB0).*im).*mask)

            # @info("Stop... high B0 corr...")
            # @infiltrate
        end

        kdata[1] = dropdims(reshape(tmp,(:,numCha,1)),dims=3)

        # ##### Temp, removing some samples... To test SNAIL correction..
        # @info("Stop... Remove some samples..")
        # @infiltrate

        # # Spliting by a fraction of the readout samples, SNAIL
        # # for a fraction of the kspace
        # increment = 1
        # split_fraction = Int(5)
        # sz = Int(ceil(params_pulseq["gen"]["ro_samples"]/split_fraction))  
        # # # for a specific number of samples
        # # split_fraction = 1
        # # sz = Int(9920)  
        # # increment = Int(params_pulseq["gen"]["ro_samples"]-sz)

        # idx = []
        # j=1

        # for i=1:Int(params_pulseq["gen"]["n"][3]* params_pulseq["spi"]["interl"])
        #     if i==1
        #         idx =[idx;(sz*(j-1))+1:(sz*(j-1)+sz)]
        #     else
        #         idx =[idx;((sz+increment)*(j-1))+1:((sz+increment)*(j-1)+sz)]
        #     end
        #     j += split_fraction
        # end
        # if params_pulseq["spi"]["type"] == 1 && split_fraction != 1
        #     idx = idx.+Int(params_pulseq["gen"]["ro_samples"] - sz)
        # end


        # @infiltrate

        # ks_traj.nodes = ks_traj.nodes[:,idx]
        # ks_traj.times = ks_traj.times[idx]
        # kdata[1] = kdata[1][idx,:]

        # ## Gaussian filter of center of k-space
        # gaussian_window = gaussian(Int(size(kdata[1],1)/params_pulseq["gen"]["n"][3]/params_pulseq["spi"]["interl"]),0.1; zerophase=true)
        # gaussian_window[1:Int(end/2)] .= 0
        # gaussian_window = repeat(gaussian_window,Int(params_pulseq["gen"]["n"][3]*params_pulseq["spi"]["interl"]))

        # kdata[1] = kdata[1].*gaussian_window

        # ###############

        # @info("Stop before recon ....")
        # @infiltrate

        if params[:do_b0_corr] && params[:do_b0_corr_seg] == 1
            Ireco = CorrectOffResonanceTimeSegmented(ks_traj,kdata,OffResonanceMap,params_pulseq,params_recon)
        elseif params[:do_b0_corr] && params[:do_b0_corr_seg] == 2
            Ireco = CorrectOffResonanceFrequencySegmented(ks_traj,kdata,OffResonanceMap,params_pulseq,params_recon)
        else
            # Create AcqData object
            acqData = AcquisitionData(ks_traj,kdata; encodingSize=(Int(params_pulseq["gen"]["ro_samples"]*params_pulseq["spi"]["interl"]),1,Int(params_pulseq["gen"]["n_ov"][3])), fov=Tuple(params_pulseq["gen"]["fov"]))
            
            # Reconstruction
            @info(string("Reconstructing Repetition # ",i_rep))
            @time Ireco = reconstruction(acqData, params_recon)
        end

        # Remove extra slices due to phase oversampling
        if params_pulseq["gen"]["n_ov"] ≠ params_pulseq["gen"]["n"]
            ph_ov_slices = params_pulseq["gen"]["n_ov"][3] - params_pulseq["gen"]["n"][3]
            ph_ov_slices = Int(ph_ov_slices/2)
            Ireco = Ireco[:,:,ph_ov_slices+1:end-ph_ov_slices]
        end

        # Save magnitude
        if params[:do_pi_recon]
            Ireco_mag = NIVolume(abs.(Ireco[:,:,:,1,1,:]))
        else
            Ireco_mag = sqrt.(sum((abs.(Ireco).^2),dims=5))
            Ireco_mag = NIVolume(abs.(Ireco_mag[:,:,:,1,:,:]))
        end

        # Save phase
        if params[:save_ph]
            Ireco_ph = NIVolume(angle.(Ireco[:,:,:,1,1,:]))
            niwrite(string(params[:path],"/recon/3d/",params[:scan],"_rep_",i_rep,save_suffix,"_ph.nii"),Ireco_ph)
        end

        niwrite(string(params[:path],"/recon/3d/",params[:scan],"_rep_",i_rep,save_suffix,".nii"),Ireco_mag)

        # Trying to clear some memory
        ccall(:malloc_trim, Cvoid, (Cint,), 0) 
        GC.gc()

    end
    
    # Merge all repetitions into one file
    for i_contrasts = 1:length(params[:contrasts])
        if params_pulseq["gen"]["seq"] == 1
            tmp_rep_recon = params[:rep_recon][1]-1+i_contrasts:2:params[:rep_recon][end]#+i_contrasts
        else
            tmp_rep_recon = params[:rep_recon]
        end 
        TimeSeries = MergeReconstrucion(params[:path],params[:scan],tmp_rep_recon,params[:do_pi_recon],params[:do_b0_corr],params[:do_coco_corr],params[:do_k0_corr],params[:do_rDORK_corr],params[:do_high_b0_corr]; TrajectoryType=params[:traj_type])
        if params[:rep_recon] == 0
            file = string(params[:path],"/recon/",params[:scan],"_",params[:contrasts][i_contrasts],save_suffix,".nii")
        else
            file = string(params[:path],"/recon/",params[:scan],"_",params[:contrasts][i_contrasts],"_r",params[:rep_recon][1],"_",params[:rep_recon][end],save_suffix,".nii")
        end
        # Creating NIfTI file, adding pixel size and TR
        TimeSeries_nii = NIVolume(TimeSeries; voxel_size=Tuple(params_pulseq["gen"]["res"].*1e3),time_step=params_pulseq["gen"]["volTR"])
        if isfile(file)
            @info ("TimeSeries exists... Re-write? y/n")
            input = readline()
            if input == "y"
                niwrite(file,TimeSeries_nii) 
            end
        else
            niwrite(file,TimeSeries_nii)
        end
    end

end
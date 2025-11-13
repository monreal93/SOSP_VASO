using Pkg

cd("/neurodesktop-storage/5T4/Alejandro/sosp_vaso")
# cd("/usr/share/5T3/Alejandro/sosp_vaso/")
Pkg.activate("./recon/")

using Revise
using Infiltrator
using MRIReco, MRIFieldmaps, MRIFiles, MRICoilSensitivities, RegularizedLeastSquares
using MAT, NIfTI , JLD, MriResearchTools
using StatsBase: mean
using Unitful: s, ustrip
using MIRTjim: jim, prompt; jim(:prompt, true)
using ImageTransformations, Interpolations
using NFFT
using FLoops
using Statistics
using DSP
using ImageMorphology
using PaddedViews
using JuMP
import Ipopt, HiGHS
# using Plots
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
include("./functions/fn_CorrectDynOffResonanceFID.jl")
include("./functions/fn_CorrectMotionFID.jl")
include("./functions/fn_CorrectOffResonance.jl")

params = Dict{Symbol, Any}()
params[:do_pi_recon] = true             # Perform PI reconstruction or direct recon
params[:do_b0_corr] = true
params[:do_b0_corr_seg] = 2   # Perform custom  0=Gridding(MRIReco), 1=Time Segmented, 2=Frequency segmented,
params[:do_rDORK_corr] = false      
params[:do_iDORK_corr] = false            # Perform Interleaf DORK (WIP)
params[:do_pDORK_corr] = false            # Perform Partition DORK (WIP)
params[:do_coco_corr] = false           # Perform concomitant field correction (WIP)  
params[:do_k0_corr] = false             # Perform K0 demodulation of raw data, only works with sk/girf trajectory (WIP)
params[:do_dyn_b0_corr] = false         # Higher order correction, FID navigators (WIP) 
params[:do_motion_corr] = false         # Motion correction, FID navigators (WIP)
params[:do_t2s_corr] = false
params[:is2d] = false
params[:isPhantom] = false              # If true, it will use params[:scan_b0] for recon
params[:rep_recon] =   4 # 169:(168*2) # 155:(154*2)            # Range of rep to recon, 0 for all rep, single number for 1 rep, or range (2:10) 
params[:traj_type] = "nom"                 # Trajectory type nominal="nom",skope="sk",poet = "poet", girf = "girf"
params[:save_ph] = false                       # Save phase of recon as nifti
params[:mcorr] = ""           # Motion correction with navigators "_mCorr"
params[:recon_order] =  1                  # Higher order recon (2,3)
params[:coil_compression] = true
b0_fq_shift = 0;                     # B0 frequency shift for B0 map...(-400 for 9.4T), +1600 for 11.7T

# Some parameters
params[:scan] = "sb_04_DS_SO_08mm_fb_rz33"            # For now: if multipe echos, include _e1.. _e2..
params[:scan_b0] = "..."           # Name of the ME-GRE to use for CS and B0map, ONLY USED WHEN SEPARATE FIELDMAP SCAN...
params[:directory] = "06242025_sb_7T"        # directory where the data is stored

# Drive where data is stored... (5T4/5T3)
path_tmp = "/neurodesktop-storage/5T4/Alejandro/sosp_vaso/"
params[:path] = string(path_tmp,"data/",params[:directory])

# Read pulseq parameters
params_pulseq = matread(string(params[:path],"/acq/",params[:scan],"_params.mat")); params_pulseq = params_pulseq["params"]

# Some constraints
if params[:recon_order] > 1; params[:traj_type] = "sk"; end
# If GIRF.. update n,n_ov and res
if params[:traj_type] == "girf"
    params_pulseq["gen"]["res"] = deepcopy(params_pulseq["girf"]["res"])
    params_pulseq["gen"]["n"] = deepcopy(params_pulseq["girf"]["n"])
    params_pulseq["gen"]["n_ov"] = deepcopy(params_pulseq["girf"]["n_ov"])
end
if params_pulseq["spi"]["type"] == 1
    params[:do_iDORK_corr] = params[:do_pDORK_corr] = false 
end
if params_pulseq["spi"]["interl"] == 1
    params[:do_iDORK_corr] = false
end

# File names
if params_pulseq["gen"]["me_gre"] > 0 && params[:isPhantom] == false
    idx = findall.("_",params[:scan])[1][1]+1:findall.("_",params[:scan])[2][1]-1
    file_name = string(string(params[:scan][1],params[:scan][idx]),"_",Int(params_pulseq["gen"]["n_ov"][1]),"_",Int(params_pulseq["gen"]["n_ov"][2]),"_",Int(params_pulseq["gen"]["n_ov"][3]),"_",params[:traj_type])
else
    idx = findall.("_",params[:scan_b0])[1][1]+1:findall.("_",params[:scan_b0])[2][1]-1
    file_name = string(string(params[:scan_b0][1],params[:scan_b0][idx]),"_",Int(params_pulseq["gen"]["n_ov"][1]),"_",Int(params_pulseq["gen"]["n_ov"][2]),"_",Int(params_pulseq["gen"]["n_ov"][3]),"_",params[:traj_type])
end
save_suffix = string("_",params[:traj_type]) 

# read spiral raw MRD file, add partition labels, format into expected array shape
file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/",params[:scan],".h5"))
@info("Reading raw data...")
rawData = RawAcquisitionData(file)

@info("Stop... save only one volume...")
@infiltrate

# Splitting raw data into fieldmap and fMRI, and calculate Drift
RepRef = 1         # DORK/OffresonanceDrift reference repetition
B0RepRef = 1       # If multiple Repetitions acquired of ME-GRE, choose wihich one to use for Drift and Sensitivity/Fieldmap
if params_pulseq["gen"]["me_gre"] == 1
    # Read pulseq parameters    
    if params[:isPhantom]
        params_me_gre_pulseq = matread(string(params[:path],"/acq/",params[:scan_b0],"_params_me_gre.mat")); params_me_gre_pulseq = params_me_gre_pulseq["params_me_gre"]
    else
        params_me_gre_pulseq = matread(string(params[:path],"/acq/",params[:scan],"_params_me_gre.mat")); params_me_gre_pulseq = params_me_gre_pulseq["params_me_gre"]
    end

    if params[:traj_type] == "girf"
        params_me_gre_pulseq["gen"]["res"] = deepcopy(params_me_gre_pulseq["girf"]["res"])
        params_me_gre_pulseq["gen"]["n"] = deepcopy(params_me_gre_pulseq["girf"]["n"])
        params_me_gre_pulseq["gen"]["n_ov"] = deepcopy(params_me_gre_pulseq["girf"]["n_ov"])
    end

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
        nav_m = mean(DSP.unwrap(angle.(nav); dims=1))
        nav_b0_m = mean(DSP.unwrap(angle.(nav_b0); dims=1))

        OffResonanceDrift = (nav_m-nav_b0_m)./params_me_gre_pulseq["gen"]["te"][1]
    end

elseif params_pulseq["gen"]["me_gre"] == 2
    if params[:isPhantom]
        file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/",params[:scan_b0],"_me_gre.h5"))
    else
        file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/",params[:scan],"_me_gre.h5"))
    end
    rawData_b0 = RawAcquisitionData(file)

    # Read pulseq parameters
    params_me_gre_pulseq = matread(string(params[:path],"/acq/",params[:scan],"_params_me_gre.mat")); params_me_gre_pulseq = params_me_gre_pulseq["params_me_gre"]

    if params[:traj_type] == "girf"
        params_me_gre_pulseq["gen"]["res"] = deepcopy(params_me_gre_pulseq["girf"]["res"])
        params_me_gre_pulseq["gen"]["n"] = deepcopy(params_me_gre_pulseq["girf"]["n"])
        params_me_gre_pulseq["gen"]["n_ov"] = deepcopy(params_me_gre_pulseq["girf"]["n_ov"])
    end

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
        nav_m = mean(DSP.unwrap(angle.(nav); dims=1))
        nav_b0_m = mean(DSP.unwrap(angle.(nav_b0); dims=1))

        # @info("Stop... B0 drift...")
        # @infiltrate

        OffResonanceDrift = (nav_m-nav_b0_m)./params_me_gre_pulseq["gen"]["te"][1]
    end
else
    # params_me_gre_pulseq = matread(string(params[:path],"/acq/",params[:scan_b0],"_params_me_gre.mat")); params_me_gre_pulseq = params_me_gre_pulseq["params_me_gre"]
end

if isfile(string(params[:path],"/acq/me_",file_name,".jld"))
    @info ("ME GRE exists... Reconstruct again? y/n")
    input = readline()
else
    input = ""
end
if input == "n"
    @info ("Loading ME GRE ...")
    recon_b0 = load(string(params[:path],"/acq/me_",file_name,".jld"))
    recon_b0 = recon_b0["recon_b0"]
    if params_pulseq["gen"]["me_gre"] != 0
        params[:b0_TE] = vec(params_me_gre_pulseq["gen"]["te"].*1e3)
    end
    recon_b0_resize = params_pulseq["gen"]["n"]
else
    @info ("Reconstructing ME GRE ...")
    # Reconstruction of ME-GRE scan (spiral from pulseq or cartesian separate scan)
    if params_pulseq["gen"]["me_gre"] > 0
        ks_traj_me_gre = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_me_gre_",params[:traj_type],".mat")); ks_traj_me_gre = ks_traj_me_gre["ks_traj_me_gre"]
        # Normalizing
        # Not sure if I want to do -1 for kx
        ks_traj_me_gre["kx"] = ks_traj_me_gre["kx"]./(maximum([abs(minimum(ks_traj_me_gre["kx"])),abs(maximum(ks_traj_me_gre["kx"]))])*2)
        ks_traj_me_gre["ky"] = ks_traj_me_gre["ky"]./(maximum([abs(minimum(ks_traj_me_gre["ky"])),abs(maximum(ks_traj_me_gre["ky"]))])*2)#.*(-1)
        ks_traj_me_gre["kz"] = ks_traj_me_gre["kz"]./(maximum([abs(minimum(ks_traj_me_gre["kz"])),abs(maximum(ks_traj_me_gre["kz"]))])*2)

        if params[:do_k0_corr] && (params[:traj_type] == "girf" || params[:traj_type] == "sk")
            k0_meas = Float64.(ks_traj_me_gre["k0"])
            k0_meas = k0_meas#.*(-1) # Do I need to do -1, I do it for ks_traj..?
            if  params[:traj_type] == "sk"
                k0_sim = matread(string(params[:path],"/acq/",params[:scan],"_k0_sim.mat")); k0_sim = k0_sim["k0"]
            end
        end
        
        ks_traj_me_gre = hcat(ks_traj_me_gre["kx"][:],ks_traj_me_gre["ky"][:],ks_traj_me_gre["kz"][:])
        ks_traj_me_gre = permutedims(ks_traj_me_gre,[2,1])

        ks_traj_me_gre = convert(AbstractMatrix{Float32},ks_traj_me_gre)
        times_me_gre = repeat(params_me_gre_pulseq["gen"]["t_vector"][:],Int(params_me_gre_pulseq["gen"]["n_ov"][3]))
        times_me_gre = convert(Vector{Float32},times_me_gre)


        ks_traj_me_gre = Trajectory(ks_traj_me_gre,1,Int(params_me_gre_pulseq["gen"]["ro_samples"]*params_me_gre_pulseq["spi"]["interl"]); 
                times=times_me_gre,TE=params_me_gre_pulseq["gen"]["TE"],AQ=params_me_gre_pulseq["gen"]["ro_time"], 
                numSlices=Int(params_me_gre_pulseq["gen"]["n_ov"][3]),circular=false)

        # Reconstruct echo by echo
        kdata = Array{Array{ComplexF32,2},3}(undef,1,1,1)
        recon_b0 = Array{ComplexF32,5}(undef,(Int.(params_me_gre_pulseq["gen"]["n_ov"])...,size(rawData.profiles[1].data,2),Int.(params_me_gre_pulseq["gen"]["me_gre_echos"])))
        for i_ech=1:params_me_gre_pulseq["gen"]["me_gre_echos"]

            b0_ech_range = Int(i_ech):Int(params_me_gre_pulseq["gen"]["me_gre_echos"]):MeGreProfiles
            b0_ech_range = b0_ech_range .+ (MeGreProfiles*(B0RepRef-1))

            rawData_b0_ech = RawAcquisitionData(rawData_b0.params,rawData_b0.profiles[b0_ech_range])

            # rawData_b0 = RawAcquisitionData(rawData.params,rawData.profiles[Int(i_ech):Int(params_me_gre_pulseq["gen"]["me_gre_echos"]):MeGreProfiles])
            rawData_b0_ech.params["TE"] = vec(params_me_gre_pulseq["gen"]["te"].*1e3)

            numRead = size(rawData_b0_ech.profiles[3].data,1) # Readout from raw data, taking 3rd in case FID
            numInterl = Int(params_me_gre_pulseq["spi"]["interl"])
            numPar = Int(params_me_gre_pulseq["gen"]["n_ov"][3]./params_me_gre_pulseq["gen"]["kz"])
            numSet = maximum(unique([rawData_b0_ech.profiles[l].head.idx.set+1 for l=1:length(rawData_b0_ech.profiles)]))
            numCha = size(rawData_b0_ech.profiles[1].data,2)
            # numRep = rawData_b0_ech.params["userParameters"]["sWipMemBlock.alFree[7]"]  # I think this is repetitions in special card
            numRep = 1
            # numEch = Int(params_me_gre_pulseq["gen"]["me_gre_echos"])
            numEch = 1

            # Adding some extra parameters
            params[:numRead] = numRead*numSet
            params[:numInterl] = numInterl
            params[:numSet] = numSet
            params[:numPar] = numPar
            params[:numCha] = numCha
            params[:numRep] = numRep
            params[:numEch] = numEch

            tmp, ndata = FormatRawData(rawData_b0_ech,params;single_rep=true,rep_format=1,MultiEchoGRE=true)

            # K0 correction
            if params[:do_k0_corr]  && params[:traj_type] == "girf"
                tmp = Correctk0(tmp,k0_meas,nothing,params)
            elseif params[:do_k0_corr]  && params[:traj_type] == "sk"
                tmp = Correctk0(tmp,k0_meas,k0_sim,params)
            end

            kdata[1] = dropdims(reshape(tmp,(:,numCha,1)),dims=3)
            
            acqData_b0 = AcquisitionData(ks_traj_me_gre,kdata; encodingSize=(Int(params_me_gre_pulseq["gen"]["ro_samples"]*params_me_gre_pulseq["spi"]["interl"]),1,Int(params_me_gre_pulseq["gen"]["n_ov"][3])), fov=Tuple(params_me_gre_pulseq["gen"]["fov"]))

            recon_b0[:,:,:,:,Int(i_ech)] = ReconSpiralMEGRE(acqData_b0,params_me_gre_pulseq)

        end

        # Remove extra slices due to phase oversampling
        if params_me_gre_pulseq["gen"]["n_ov"] ≠ params_me_gre_pulseq["gen"]["n"]
            ph_ov_slices = (params_me_gre_pulseq["gen"]["n_ov"][3]*params_me_gre_pulseq["gen"]["kz"]) - params_me_gre_pulseq["gen"]["n"][3]
            ph_ov_slices = Int(ph_ov_slices/2)
            recon_b0 = recon_b0[:,:,ph_ov_slices+1:end-ph_ov_slices,:,:,:]
        end
                
        params[:b0_TE] = vec(params_me_gre_pulseq["gen"]["te"].*1e3)
    else
        # read ME-GRE MRD data and recon
        file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/",params[:scan_b0],".h5"))
        acqData_b0 = AcquisitionData(file)
        rawData_b0 = RawAcquisitionData(file)

        if params_pulseq["gen"]["field_strength"] == 7im
            recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=false)
        elseif params_pulseq["gen"]["field_strength"] == 7
            recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=true)
        elseif params_pulseq["gen"]["field_strength"] == 9
            recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=true) 
        elseif params_pulseq["gen"]["field_strength"] == 11
            recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=true)
        end
    end

    # # Remove extra slices due to phase oversampling
    # if params_pulseq["gen"]["me_gre"] > 0
    #     if (params_me_gre_pulseq["gen"]["n_ov"] ≠ params_me_gre_pulseq["gen"]["n"])
    #         ph_ov_slices = (params_me_gre_pulseq["gen"]["n_ov"][3]*params_me_gre_pulseq["gen"]["kz"]) - params_me_gre_pulseq["gen"]["n"][3]
    #         ph_ov_slices = Int(ph_ov_slices/2)
    #         recon_b0 = recon_b0[:,:,ph_ov_slices+1:end-ph_ov_slices,:,:]
    #     end
    # end

    save(string(params[:path],"/acq/me_",file_name,".jld"),"recon_b0",recon_b0)

    if params_pulseq["gen"]["me_gre"] == 0
        # # AMM: For now I add them manually... Fix it..
        @info("Manually write the TE for B0...")
        @infiltrate
        params[:b0_TE] = vec([2.5,4.5,7])
    end

    # Saving ME GRE for reference in tmp folder
    b0_nii = recon_b0[:,:,:,:,:]
    b0_nii = dropdims(sqrt.(sum(abs.(b0_nii).^2, dims=4)), dims=4)
    # Note AMM (03062025): seems like I need to reshape to 60 in z... also need to circshift 
    # circshift(b0_nii,(0,2,-2,0));
    @info("stop... recon_b0...")
    @infiltrate
    if params_pulseq["gen"]["me_gre"] > 0
        if params_pulseq["gen"]["n_ov"][3] != params_me_gre_pulseq["gen"]["n_ov"][3]
            recon_b0_resize = params_me_gre_pulseq["gen"]["fov"] ./ params_pulseq["gen"]["res"]
            recon_b0_resize[1:2] = round.(params_me_gre_pulseq["gen"]["fov"][1:2] ./ params_pulseq["gen"]["res"][1:2])  # oP1
            recon_b0_resize[3] = round.(params_me_gre_pulseq["gen"]["fov"][3] ./ params_pulseq["gen"]["res"][3]./2)*2 
            # recon_b0_resize[1] = 258;
            # recon_b0_resize = [params_pulseq["gen"]["n"][1:2]... round.(params_me_gre_pulseq["gen"]["fov"][3] ./ params_pulseq["gen"]["res"][3])] # oP 2
            b0_nii = imresize(b0_nii,Int(recon_b0_resize[1]),Int(recon_b0_resize[2]),Int(recon_b0_resize[3]))
            crop_size = Int(ceil((recon_b0_resize[3]-params_pulseq["gen"]["n"][3])/2))
            if crop_size > 0
                b0_nii = b0_nii[:,:,crop_size+1:end-crop_size+1,:]
            end
        else
            recon_b0_resize = params_pulseq["gen"]["n_ov"]
            b0_nii = imresize(b0_nii,Int(recon_b0_resize[1]),Int(recon_b0_resize[2]),Int(recon_b0_resize[3])) 
        end
    else
        recon_b0_resize = params_pulseq["gen"]["n_ov"]
        b0_nii = imresize(b0_nii,Int(recon_b0_resize[1]),Int(recon_b0_resize[2]),Int(recon_b0_resize[3]))
    end

    # b0_nii = imresize(b0_nii,Int(params_pulseq["gen"]["n_ov"][1]),Int(params_pulseq["gen"]["n_ov"][2]),Int(params_pulseq["gen"]["n_ov"][3]))
    # if rawData.profiles[1].head.read_dir[1] == 1 && rawData.profiles[1].head.phase_dir[3] == 1 && rawData.profiles[1].head.slice_dir[2] == 1     # Coronal
    #     quatern_b=0.7; quatern_c=0; quatern_d=0;
    # elseif rawData.profiles[1].head.read_dir[1] == 1 && rawData.profiles[1].head.phase_dir[3] == 1 && rawData.profiles[1].head.slice_dir[2] == 1
    #     # ToDo
    # elseif rawData.profiles[1].head.read_dir[1] == 1 && rawData.profiles[1].head.phase_dir[3] == 1 && rawData.profiles[1].head.slice_dir[2] == 1
    #     # ToDo
    # end
    # b0_nii = NIVolume(b0_nii; voxel_size=Tuple(params_pulseq["gen"]["res"].*1e3), 
        # qfac=Float32(-1), quatern_b=quatern_b ,quatern_c=quatern_c ,quatern_d=quatern_d ,
        # qoffset_x = rawData.profiles[1].head.position[1]/2, qoffset_y=rawData.profiles[1].head.position[2]/2, qoffset_z=rawData.profiles[1].head.position[3]/2,
        # )
    # orientation=Float32.([params_pulseq["gen"]["res"][1].*1e3 0 0 0;0 0 params_pulseq["gen"]["res"][3].*1e3 0; 0 params_pulseq["gen"]["res"][2].*1e3 0 0])
    # b0_nii = NIVolume(b0_nii; voxel_size=Tuple(params_pulseq["gen"]["res"].*1e3))

    b0_nii = NIVolume(b0_nii)
    niwrite(string(params[:path],"/tmp/",params[:scan],"_",params[:traj_type],".nii"),b0_nii)
end

recon_b0_sos = dropdims(sqrt.(((sum(abs.(recon_b0).^2; dims=4)))), dims=4)

# Getting some parameters from the raw data
numRead = size(rawData.profiles[3].data,1) # Readout from raw data, taking 3rd in case FID
# numRead = Int(params_pulseq["gen"]["ro_samples"])  # Readout from pulseq params
numInterl = Int(params_pulseq["spi"]["interl"])
# numPar = maximum(unique([rawData.profiles[l].head.idx.kspace_encode_step_2+1 for l=1:length(rawData.profiles)]))
numPar = Int(params_pulseq["gen"]["n_ov"][3])
numSet = maximum(unique([rawData.profiles[l].head.idx.set+1 for l=1:length(rawData.profiles)]))
numCha = size(rawData.profiles[1].data,2)
# numRep = rawData.params["userParameters"]["sWipMemBlock.alFree[7]"]  # I think this is repetitions in special card
numRep = 4
numEch = Int(params_pulseq["gen"]["echos"])
# numRep = Int(length(rawData.profiles)/numPar/numSet)
if params_pulseq["gen"]["seq"] == 3
    numRep = numRep * Int(length(params_pulseq["gen"]["multi_te"]))
end

# Adding some extra parameters
params[:numRead] = numRead*numSet
params[:numInterl] = numInterl
params[:numSet] = numSet
params[:numPar] = numPar
params[:numCha] = numCha
params[:numRep] = numRep
params[:numEch] = numEch
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
elseif params_pulseq["gen"]["seq"] == 3
    params[:contrasts] = ["b1"]
    for i_contr in 1:length(params_pulseq["gen"]["multi_te"]);params[:contrasts]=cat(params[:contrasts],string("b",i_contr),dims=1); end
elseif params_pulseq["gen"]["seq"] == 4
    # BOLD
    params[:contrasts] = ["b"]
end

# AMM: Temp... trying to mask CS
@info("stop.. temp...")
@infiltrate
# mask = recon_b0_sos[:,:,:,1]
# mask[abs.(mask) .< 0.5e-7] .= 0
# mask[abs.(mask) .!= 0] .= 1
# mask = isone.(mask)
# recon_b0 .*= mask

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
        # SensitivityMap = CalculateSensitivityMap(recon_b0,Tuple(Int.(params_pulseq["gen"]["n_ov"].*[1 1 params_pulseq["gen"]["kz"]]))) # Original
        SensitivityMap = CalculateSensitivityMap(recon_b0,Tuple(Int.(recon_b0_resize)))
        save(string(params[:path],"/acq/cs_",file_name,".jld"),"SensitivityMap",SensitivityMap)
    end
end

# Calculate B0 map
if params[:do_b0_corr] || (params[:traj_type] == "sk" && params[:recon_order] > 1)
    if params[:do_dyn_b0_corr] && params[:do_b0_corr_seg] == 0
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
        global OffResonanceMap = OffResonanceMap["OffResonanceMap"]
    else
        @info("Calculating Offresonance Maps....")
        global OffResonanceMap = CalculateOffresonanceMap(recon_b0,SensitivityMap,params[:b0_TE])
        save(string(params[:path],"/acq/b0_",file_name,".jld"),"OffResonanceMap",OffResonanceMap)
    end
end

# Calculate concomitant field
if params[:do_coco_corr] && params[:do_b0_corr]
        save_suffix = string(save_suffix,"_co")
        if isfile(string(params[:path],"/acq/co_",file_name,".jld"))
            @info ("Concomitant map exists... Re-calculate? y/n")
            input = readline()
        else
            input = ""
        end
        if input == "n"
            @info ("Loading Concomitant maps ...")
            CocoFieldMap = load(string(params[:path],"/acq/co_",file_name,".jld"))
            CocoFieldMap = CocoFieldMap["CocoFieldMap"]
            OffResonanceMap .= OffResonanceMap .+ (CocoFieldMap.*im)
        else
            @info("Calculating Concomitant Maps....")
            # phase,read,slice
            RotMat =vcat(transpose(collect(rawData_b0.profiles[1].head.read_dir)),
                    transpose(collect(rawData_b0.profiles[1].head.phase_dir)),
                    transpose(collect(rawData_b0.profiles[1].head.slice_dir)))
            CenterPosition = (rawData_b0.profiles[end].head.position .+ rawData_b0.profiles[1].head.position)./2 .* 1e-3
            CocoFieldMap = CalculateConcomitantFieldMap(RotMat,CenterPosition, params_pulseq)
            # CocoFieldMap = CocoFieldMap .* (2*π)
            # CocoFieldMap = CocoFieldMap .- Float32(10*2π) # Ading some extra Hz, maybe due to drift?
            # I need to reverse dims=3,2,1, maybe rot matrix is wrong?
            CocoFieldMap = reverse(CocoFieldMap,dims=(2)) # maybe just for 7T?
            # OffResonanceMap .= OffResonanceMap .* (π/2)
            # Sometimes - works better, sometimes +
            OffResonanceMap .= OffResonanceMap .- (CocoFieldMap.*im)
            # OffResonanceMap .= OffResonanceMap .- (14*2*π*im)
            mask = SensitivityMap[:,:,:,1]
            mask[abs.(mask) .> 0] .= 1
            mask = isone.(mask)
            OffResonanceMap .= OffResonanceMap .* mask
            save(string(params[:path],"/acq/co_",file_name,".jld"),"CocoFieldMap",CocoFieldMap)
        end
end

# Creating Brain mask and Brain+skull mask
# AMM: Ideally I would like to do this either when I recon ME-GRE or calculate Sensitivity/OffResonance
# cs = imresize(SensitivityMap,size(recon_b0)[1:4])
# recon_b0_cc = sum(conj.(cs) .* recon_b0; dims=4)
# recon_b0_cc = dropdims(recon_b0_cc, dims=4)
# phase=angle.(recon_b0_cc[:,:,:,1])
# mag=abs.(recon_b0_cc[:,:,:,1])
# brain = brain_mask(robustmask(romeovoxelquality(phase; mag); threshold=0.6))
# brain_skull_mask = Float32.(dilate(brain; r=15))
# brain_skull_mask = imresize(brain_skull_mask, size(SensitivityMap)[1:3])
# if params[:do_pi_recon]; SensitivityMap = SensitivityMap .* brain_skull_mask; end
# if params[:do_b0_corr]; OffResonanceMap = OffResonanceMap .* brain_skull_mask; end

# @info("Stop... temp...")
# @infiltrate

# Normalize k-space trajectory and create Trajectory object
if params[:traj_type] == "sk"
    ks_traj = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_",params[:traj_type],".mat")); ks_traj = ks_traj["ks_traj"]
    # ToDo: Sync Skope trajectory to raw data...

    ## Approach 2) same structure as nominal
    # Normalizing
    ks_traj["kx"] = ks_traj["kx"]./(maximum([abs(minimum(ks_traj["kx"])),abs(maximum(ks_traj["kx"]))])*2)
    ks_traj["ky"] = ks_traj["ky"]./(maximum([abs(minimum(ks_traj["ky"])),abs(maximum(ks_traj["ky"]))])*2)#.*(-1)
    ks_traj["kz"] = ks_traj["kz"]./(maximum([abs(minimum(ks_traj["kz"])),abs(maximum(ks_traj["kz"]))])*2)
    # k0 = ks_traj["k0"][:]
    ks_traj = hcat(ks_traj["kx"][:],ks_traj["ky"][:],ks_traj["kz"][:])
    ks_traj = permutedims(ks_traj,[2,1])

else
    ks_traj = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_",params[:traj_type],".mat")); ks_traj = ks_traj["ks_traj"]
    # Normalizing
    ks_traj["kx"] = ks_traj["kx"]./(maximum([abs(minimum(ks_traj["kx"])),abs(maximum(ks_traj["kx"]))])*2)
    ks_traj["ky"] = ks_traj["ky"]./(maximum([abs(minimum(ks_traj["ky"])),abs(maximum(ks_traj["ky"]))])*2)#.*(-1)
    ks_traj["kz"] = ks_traj["kz"]./(maximum([abs(minimum(ks_traj["kz"])),abs(maximum(ks_traj["kz"]))])*2)
    ks_traj = hcat(ks_traj["kx"][:],ks_traj["ky"][:],ks_traj["kz"][:])
    ks_traj = permutedims(ks_traj,[2,1])
end

# K0 from skope
if params[:do_k0_corr] && params[:traj_type] == "sk"
    k0_meas = matread(string(params[:path],"/acq/",params[:scan],"_k0_sk.mat")); k0_meas = k0_meas["k0"]# [:]
    k0_sim = matread(string(params[:path],"/acq/",params[:scan],"_k0_sim.mat")); k0_sim = k0_sim["k0"]
elseif params[:do_k0_corr] && params[:traj_type] == "girf"
    @info("Stop... K0..")
    @infiltrate
    k0_meas = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_",params[:traj_type],".mat")); k0_meas = k0_meas["ks_traj"]
    k0_meas = Float64.(k0_meas["k0"])
    k0_meas = k0_meas#.*(-1) # Do I need to do -1, I do it for ks_traj..?
end

ks_traj = convert(AbstractMatrix{Float32},ks_traj)
times = repeat(params_pulseq["gen"]["t_vector"][:],numPar)
times = convert(Vector{Float32},times)

ks_traj = Trajectory(ks_traj,1,Int(params_pulseq["gen"]["ro_samples"]*params_pulseq["spi"]["interl"]); 
        times=times,TE=params_pulseq["gen"]["TE"],AQ=params_pulseq["gen"]["ro_time"], 
        numSlices=Int(params_pulseq["gen"]["n_ov"][3]),circular=false)


# Setting repetition label
rawData = SetRepetitionLabel(rawData,params)

# K0 correction
if params[:do_k0_corr]
    save_suffix = string(save_suffix,"_k0")
end

# Cropping Sensitivity and OffResonance maps
if params_pulseq["gen"]["me_gre"] > 0
    # If me_gre FOV > FOV
    if params[:do_pi_recon] && (params_pulseq["gen"]["fov_ov"][3] != params_me_gre_pulseq["gen"]["fov_ov"][3])
        # crop_size = Int(ceil((size(SensitivityMap,3) - (params_pulseq["gen"]["n"][3]))/2))
        # SensitivityMap_crop = SensitivityMap[:,:,crop_size+1:end-crop_size+1,:]
        crop_size = Int.(floor.((collect(size(SensitivityMap)[1:3]) .- Tuple(Int.(params_pulseq["gen"]["n"])))))
        crop_size1 = Int.(crop_size.-floor.(crop_size./2))
        crop_size2 = Int.(crop_size.-crop_size1)
        SensitivityMap_crop = SensitivityMap[crop_size1[1]+1:end-crop_size2[1],crop_size1[2]+1:end-crop_size2[2],crop_size1[3]+1:end-crop_size2[3],:]
        if params[:do_b0_corr] 
            # OffResonanceMap_crop = OffResonanceMap[:,:,crop_size+1:end-crop_size+1]
            OffResonanceMap_crop = OffResonanceMap[crop_size1[1]+1:end-crop_size2[1],crop_size1[2]+1:end-crop_size2[2],crop_size1[3]+1:end-crop_size2[3]]
        end
    elseif params[:do_pi_recon] && (params_pulseq["gen"]["fov_ov"][3] == params_me_gre_pulseq["gen"]["fov_ov"][3])
        SensitivityMap_crop = SensitivityMap
        if params[:do_b0_corr] 
            OffResonanceMap_crop = OffResonanceMap
        end
    end
else
    SensitivityMap_crop = SensitivityMap
end

# Reconstruction parameters
params_recon = Dict{Symbol, Any}()
params_recon = merge(defaultRecoParams(), params_recon)
if params[:do_pi_recon]
    params_recon[:reco] = "multiCoil"
    params_recon[:senseMaps] = SensitivityMap_crop
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
                # OffResonanceMap = OffResonanceMap .+ ComplexF32(OffResonanceDrift.*im)
                # OffResonanceMap = OffResonanceMap .- ComplexF32(300*im) # Maybe I am missing some Hz/rad..
                # OffResonanceMap = OffResonanceMap .- ComplexF32(maximum(imag.(OffResonanceMap)).*im)  # AMM: Should I have positive values?? 
            else
                # OffResonanceMap = OffResonanceMap .- ComplexF32(abs(OffResonanceDrift).*im)
            end
        end
    end
    ######## ******* AMM: Temp: Adding some rad to fieldmap.. test...
    # OffResonanceMap[imag.(OffResonanceMap).<-1000] .= 0
    # OffResonanceMap[imag.(OffResonanceMap).>1000] .= 0
    OffResonanceMap_crop = OffResonanceMap_crop .+ ComplexF32(b0_fq_shift*im) # (-400 for 9.4T), +1600 for 11.7T
    # OffResonanceMap = OffResonanceMap.*2.5
    ######## ******* AMM: Temp: Adding some rad to fieldmap.. test...
    if params[:do_b0_corr_seg] == 0
        params_recon[:correctionMap] = OffResonanceMap_crop  # Original
    end
end
if params[:recon_order] > 1 || (params[:traj_type] == "sk" && params[:recon_order] > 1)
    params_recon[:reco] = "expandedSignalModel"
    params_recon[:higherOrderTrajectory] = convert(Matrix{Float32},ks_traj_high)
    params_recon[:senseMaps] = SensitivityMap
    params_recon[:correctionMap] = OffResonanceMap_crop
end

params_pulseq["gen"]["n_ov"] = Int.(params_pulseq["gen"]["n_ov"])
params_pulseq["gen"]["n"] = Int.(params_pulseq["gen"]["n"])
params_recon[:reconSize] = Tuple(params_pulseq["gen"]["n"])
# params_recon[:regularization] = "L2"
# params_recon[:λ] = 1e3  # 1.e-2
# params_recon[:iterations] = 40 # (10)
# params_recon[:solver] = "cgnr" # cgnr (L2), fista (L1), admm(L1)
# params_recon[:method] = "nfft" # nfft, leastsquare
### New definitions
params_recon[:solver] = ADMM
params_recon[:reg] = L1Regularization(1e-3) # (1e-3)
# params_recon[:sparseTrafo] = "Wavelet"
params_recon[:iterations] = 40

# FFT Shift
fftShift = (rawData.profiles[1].head.position .- rawData_b0.profiles[1].head.position)./Float32(2)

# Getting navigator range
if Int(params_pulseq["gen"]["fid_nav"]) == 1
    params[:fid_nav_ΔTE] = params_pulseq["gen"]["fid_nav_te"][2]-params_pulseq["gen"]["fid_nav_te"][1]
    # nav_range = 6:50 # Original
    nav_range = 40:60
else 
    if params_pulseq["spi"]["type"] == 0
        nav_range = 6:20
    elseif params_pulseq["spi"]["type"] == 1
       nav_range = Int.(params_pulseq["gen"]["ro_samples"]-19:params_pulseq["gen"]["ro_samples"]-5)
    end
end 

# repetition DORK correction
nav_ref, ndata_ref = FormatRawData(rawData,params;single_rep=true,rep_format=RepRef,fid_nav=Int(params_pulseq["gen"]["fid_nav"]))
if params[:do_rDORK_corr]
    save_suffix = string(save_suffix,"_rDORK")
    # RepRef = 3         # DORK reference repetition
    # ndata format = (numRead,numInterl,numPar,numCha,numRep,numFID)
    # nav_ref, ndata_ref = FormatRawData(rawData,params;single_rep=true,rep_format=RepRef,fid_nav=Int(params_pulseq["gen"]["fid_nav"]))

    if params_pulseq["gen"]["fid_nav"] == 1 && params_pulseq["gen"]["seq"] ≠ 1
        # Separate Navigator module
        # nav_ref = mean(ndata_ref[nav_range,1,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1,1],dims=2)  # Mean over channels
        # nav_ref = ndata_ref[nav_range,1,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1,1]   # Only center partition
        nav_ref = ndata_ref[nav_range,1,:,:,1,1]  
    elseif  params_pulseq["gen"]["fid_nav"] == 1 && params_pulseq["gen"]["seq"] == 1
        # VASO
        if params_pulseq["gen"]["seq"] == 1
            RepRef = RepRef+1         # DORK reference repetition
            nav_ref1, ndata_ref1 = FormatRawData(rawData,params;single_rep=true,rep_format=RepRef,fid_nav=Int(params_pulseq["gen"]["fid_nav"]))
            if Int(params_pulseq["gen"]["fid_nav"]) == 1
                # nav_ref = mean(deepcopy(ndata_ref[nav_range,:,:,:,1,1]),dims=4)   # Taking only the first FID, mean channels
                # nav_ref1 = mean(deepcopy(ndata_ref1[nav_range,:,:,:,1,1]),dims=4)   # Taking only the first FID, mean channels
                nav_ref = deepcopy(ndata_ref[nav_range,:,:,:,1,1])   # Taking only the first FID
                nav_ref1 = deepcopy(ndata_ref1[nav_range,:,:,:,1,1])   # Taking only the first FID
            end
            # if params[:kz_enc] == 0  # Linear
            #     nav_ref1 = nav_ref1[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
            #     nav_ref1 = mean(nav_ref1, dims=2)
            # elseif params[:kz_enc] == 1 # Center-out
            #     nav_ref1 = nav_ref1[nav_range,1,:,1]
            #     nav_ref1 = mean(nav_ref1, dims=2)
            # end
        end
    else
        # No separate Navigator
        if params[:kz_enc] == 0  # Linear
            nav_ref = nav_ref[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
            nav_ref = mean(nav_ref, dims=2)      # Mean over channels
        elseif params[:kz_enc] == 1 # Center-out
            nav_ref = nav_ref[nav_range,1,:,1]
            nav_ref = mean(nav_ref, dims=2)     # Mean over channels
        end
    end
end

# interleave DORK correction
if params[:do_iDORK_corr]
    save_suffix = string(save_suffix,"_iDORK")
end

# partition DORK correction
if params[:do_pDORK_corr]
    save_suffix = string(save_suffix,"_pDORK")
end

# Higher order B0 correction, FID navigators approach
if params[:do_dyn_b0_corr] && params[:do_b0_corr]
    # Calibration matrix
    DynOffResonanceCalibrationMatrix, sh_basis = getDynamicOffResonanceCalibrationMatrix(params,OffResonanceMap,recon_b0)
end

# Motion correction, FID navigators approach
if params[:do_motion_corr]
    save_suffix = string(save_suffix,"_mCorr")
    @info("Calculating motion calibration matrix....")
    # Calibration matrix
    nav_ref_mot = ndata_ref[nav_range,1,:,:,1,1]
    nav_ref_mot = abs.(dropdims(mean(nav_ref_mot, dims=(1,2)), dims=(1,2)))
    MotionCalibrationMatrix = getMotionCalibrationMatrix(recon_b0,SensitivityMap,nav_ref_mot; calib_steps=200)
    # MotionCalibrationMatrix = zeros(31,6) # Temp.. to test manual moiton params..
end

# Reconstruct repetition by repetition....
kdata = Array{Array{ComplexF32,2},3}(undef,1,1,1)
if params[:rep_recon] == 0:0 || params[:rep_recon] == 0 
    params[:rep_recon] = 1:params[:numRep]
end

@time begin
    for i_rep=params[:rep_recon]

        tmp, ndata = FormatRawData(rawData,params;single_rep=true,rep_format=i_rep,fid_nav=Int(params_pulseq["gen"]["fid_nav"]))

        # FFT Shift
        # ToDo: Use correct value for shift... how do I get it?
        tmp = CorrectFFTShift(tmp,fftShift,ks_traj.nodes,params)

        # K0 correction
        if params[:do_k0_corr] && params[:traj_type] == "sk"
            tmp = Correctk0(tmp,k0_meas,k0_sim,params)
        elseif params[:do_k0_corr] && params[:traj_type] == "girf"
            tmp = Correctk0(tmp,k0_meas,nothing,params)
        end

        if Int(params_pulseq["gen"]["fid_nav"]) == 1
            # nav = mean(ndata[nav_range,:,:,:,1,1],dims=4)      # Taking just first FID, mean over ch
            # nav = ndata[nav_range,:,:,:,:,1]                    # Taking just the first FID
            nav = ndata[nav_range,:,:,:,:,:]
            # nav = dropdims(nav, dims=5)
        else
            nav = mean(tmp[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1],dims=2)      # Taking just first FID, mean over ch
        end

        # repetition DORK correction
        if params[:do_rDORK_corr]
            if params_pulseq["gen"]["seq"] == 1  # VASO and BOLD separately
                if isodd(i_rep)  # VASO
                    tmp = CorrectRepetitionDORK(tmp,nav,nav_ref,params)
                else            # BOLD
                    tmp = CorrectRepetitionDORK(tmp,nav,nav_ref1,params)
                end
            else
                tmp = CorrectRepetitionDORK(tmp,nav,nav_ref,params)
            end
        end

        # interleaf DORK correction
        if params[:do_iDORK_corr]
            tmp = CorrectInterleafDORK(tmp,nav,params)
        end

        # partition DORK correction
        if params[:do_pDORK_corr] && Int(params_pulseq["gen"]["fid_nav"]) == 1
            tmp = CorrectPartitionDORK(tmp,nav,params)
        end

        # Dynamic Off-resonance correction.. FID navigators
        # AMM:Should I do it before or after DORK.. (I think before)
        if params[:do_dyn_b0_corr]
            @info("Stop... Minimization problem...")
            @infiltrate
            # if isnothing(ndata)
            #     nav = mean(tmp[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1],dims=1)      # Taking just first FID
            # else
            #     nav = mean(ndata[nav_range,:,:,1,1],dims=2)      # Taking just first FID
            # end
            # nav = mean(nav[:,1,1,:,:,1,1],dims=1)[:]
            # nav = nav'

            # Solve ominimization problem... 
            y_jn = mean(nav[:,1,1,:,:,1,1],dims=1)[:]
            s_j0 = recon_b0[:,:,:,:,1] # 1st echo
            s_j0 = reshape(s_j0,(:,numCha))
            s_j0 = permutedims(s_j0,[2,1]) # numCha x numPixels
            b_0 = ones(size(sh_basis,2)) # + im.*ones(size(sh_basis,2))

            f_bn(b_n) = s_j0 * exp.(im.*2pi.* 42.58e6.*params_pulseq["gen"]["fid_nav_te"][1].*(sh_basis * b_n)) 

            model = Model(Ipopt.Optimizer)
            set_silent(model)
            # @variable(model,y_jn)
            @variable(model,b_n[i=1:9],start =b_0[i]) # b_n is the coefficients of spherical harmonics
            @objective(model,Min, sum((y_jn - f_bn(b_n)).^2))
            # @objective(model,Min, (y_jn - f_bn(b_n)))

            # YY = vcat(real.(nav),imag.(nav))
            # AA = vcat(real.(DynOffResonanceCalibrationMatrix),imag.(DynOffResonanceCalibrationMatrix))

            # # YY = [real.(nav) imag.(nav)]
            # # AA = [real.(DynOffResonanceCalibrationMatrix) imag.(DynOffResonanceCalibrationMatrix)]
            
            # b = AA \ YY
            # δB0 = sh_basis * (b)

            # δB0 = reshape(δB0,Tuple(params_pulseq["gen"]["n_ov"]))

            # # δB0 = reverse(δB0,dims=(3))
            # mask = SensitivityMap[:,:,:,1]
            # mask[abs.(mask) .> 0] .= 1
            # mask = isone.(mask)

            # @info("Stop... high B0 corr...")
            # @infiltrate

            # OffResonanceMap .-= abs.(δB0).*im.*mask

            # if params[:do_b0_corr] && params[:do_b0_corr_seg] == 0
            #     params_recon[:correctionMap] = deepcopy(OffResonanceMap)
            # end


        end

        # Motion correction.. FID navigators
        if params[:do_motion_corr]
            # if isnothing(ndata)
            #     nav = mean(tmp[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,:,1],dims=1)      # Taking just first FID
            # else
            #     nav = mean(ndata[end-2:end,:,:,:,1,1],dims=2)      # Taking just first FID
            # end
            nav_mot = ndata[nav_range,1,:,:,1,1]
            nav_mot = abs.(dropdims(mean(nav_mot, dims=(1,2)), dims=(1,2)))
            tmp, ksTraj = applyMotionCalibrationMatrix(tmp,ks_traj.nodes,nav_mot,MotionCalibrationMatrix)
            ks_traj.nodes = ksTraj
        end

        kdata[1] = dropdims(reshape(tmp,(:,numCha,1)),dims=3)

        # Coil Compression 
        if params[:coil_compression] && params[:do_pi_recon]
            numVC = 16
            kdata[1], ccMat = softwareCoilCompression(kdata[1],numVC)
            params_recon[:senseMaps] = applyCoilCompressionSensitivityMaps(SensitivityMap_crop, ccMat, numVC)
        end

        @info("Stop before recon ....")
        @infiltrate

        # Reconstruction
        @info(string("Reconstructing Repetition # ",i_rep))
        if params[:do_b0_corr] && params[:do_b0_corr_seg] == 1
            @time Ireco = CorrectOffResonanceTimeSegmented(ks_traj,kdata,OffResonanceMap_crop,params_pulseq,params_recon)
        elseif params[:do_b0_corr] && params[:do_b0_corr_seg] == 2
            @time Ireco = CorrectOffResonanceFrequencySegmented(ks_traj,kdata,OffResonanceMap_crop,params_pulseq,params_recon)
        else
            # Create AcqData object
            acqData = AcquisitionData(ks_traj,kdata; encodingSize=(Int(params_pulseq["gen"]["ro_samples"]*params_pulseq["spi"]["interl"]),1,Int(params_pulseq["gen"]["n_ov"][3])), fov=Tuple(params_pulseq["gen"]["fov"]))
            @time Ireco = reconstruction(acqData, params_recon)
        end

        # Remove extra slices due to phase oversampling
        @info("Stop... final recon...")
        @infiltrate
        # if params_pulseq["gen"]["n_ov"] ≠ params_pulseq["gen"]["n"]
        #     ph_ov_slices = (params_pulseq["gen"]["n_ov"][3]*params_pulseq["gen"]["kz"]) - params_pulseq["gen"]["n"][3]
        #     ph_ov_slices = Int(ph_ov_slices/2)
        #     Ireco = Ireco[:,:,ph_ov_slices+1:end-ph_ov_slices,:,:,:]
        # end

        # Save magnitude
        if params[:do_pi_recon]
            Ireco_mag = NIVolume(abs.(Ireco[:,:,:,1,1,:]))
        else
            Ireco_mag = sqrt.(sum((abs.(Ireco).^2),dims=5))
            Ireco_mag = NIVolume(abs.(Ireco_mag[:,:,:,1]))
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
        elseif params_pulseq["gen"]["seq"] == 3
            tmp_rep_recon = params[:rep_recon][1]-1+i_contrasts:length(params_pulseq["gen"]["multi_te"]):params[:rep_recon][end]#+i_contrasts
        else
            tmp_rep_recon = params[:rep_recon]
        end 
        TimeSeries = MergeReconstrucion(params[:path],params[:scan],tmp_rep_recon,params[:do_pi_recon],params[:do_b0_corr],
                    params[:do_coco_corr],params[:do_k0_corr],params[:do_rDORK_corr],params[:do_dyn_b0_corr]; 
                    TrajectoryType=params[:traj_type], B0CorrectionType=params[:do_b0_corr_seg])
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

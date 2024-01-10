using Pkg

cd("/usr/share/5T4/Alejandro/sosp_vaso/")
Pkg.activate("./recon/")

using Revise
using Infiltrator
using MRIReco, MRIFiles, MRICoilSensitivities, MRIFieldmaps
using MAT, NIfTI , JLD #, MriResearchTools
using StatsBase: mean
using Unitful: s, ustrip
using MIRTjim: jim, prompt; jim(:prompt, true)
using ImageTransformations
using NFFT
using FLoops

using FFTW
FFTW.set_provider!("mkl")
# FFTW.set_provider!("fftw")

include("./functions/fn_save_nii.jl")
include("./functions/fn_Utilities.jl")
include("./functions/fn_PrepareMRD.jl")
include("./functions/fn_ReconstructCartesian.jl")
include("./functions/fn_CorrectDORK.jl")
include("./functions/fn_CalculateSensitivityOffresonanceMaps.jl")
include("./functions/fn_calculateSphericalHarmonics.jl")
include("../recon/functions/fn_motionCorrection.jl")

params = Dict{Symbol, Any}()
params[:do_pi_recon] = true             # Perform PI reconstruction or direct recon
params[:do_b0_corr] = false 
params[:do_rDORK_corr] = false          
params[:do_t2s_corr] = false
params[:is2d] = false
params[:rep_recon] = 1:1                      # Range of rep to recon, 0 for all rep, ex. (5:5)- rep 5 
params[:contrasts] = ["b"]                  # Contrasts to recon v,b,abc or all
params[:traj_type] = "nom"                 # Trajectory type nominal="nom",skope="sk",poet = "poet", corrected = "nom_corr"
params[:save_ph] = 0                       # Save phase of recon as nifti
params[:mcorr] = ""           # Motion correction with navigators "_mCorr"
params[:recon_order] =  1;                  # Higher order recon (2,3)

# Some parameters
params[:scan] = "sb_43"            # For now: if multipe echos, include _e1.. _e2..
params[:scan_b0] = "s41"           # Name of the ME-GRE to use for CS and B0map
params[:directory] = "10162023_sb_7Ti"        # directory where the data is stored

path_tmp = string(pwd(),'/')
params[:path] = string(path_tmp,"data/",params[:directory])

# Read pulseq parameters
params_pulseq = matread(string(params[:path],"/acq/",params[:scan],"_params.mat")); params_pulseq = params_pulseq["params"]

# read spiral raw MRD file, add partition labels, format into expected array shape
file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/",params[:scan],".h5"))
rawData = RawAcquisitionData(file)

# Getting some parameters from the raw data
numRead = size(rawData.profiles[1].data,1) 
numInterl = Int(params_pulseq["spi"]["interl"])
# numPar = maximum(unique([rawData.profiles[l].head.idx.kspace_encode_step_2+1 for l=1:length(rawData.profiles)]))
numPar = Int(params_pulseq["gen"]["n_ov"][3])
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
params[:acq_times] = params_pulseq["gen"]["t_vector"]
params[:reconSize] = (params_pulseq["gen"]["n_ov"][1],params_pulseq["gen"]["n_ov"][2],params_pulseq["gen"]["n_ov"][3])
params[:FieldOfView] = params_pulseq["gen"]["fov"]

# Normalize k-space trajectory and create Trajectory object
ks_traj = matread(string(params[:path],"/acq/",params[:scan],"_ks_traj_",params[:traj_type],".mat")); ks_traj = ks_traj["ks_traj"]
ks_traj["kx"] = ks_traj["kx"]./(maximum([abs(minimum(ks_traj["kx"])),abs(maximum(ks_traj["kx"]))])*2)
ks_traj["ky"] = ks_traj["ky"]./(maximum([abs(minimum(ks_traj["ky"])),abs(maximum(ks_traj["ky"]))])*2).*(-1)
ks_traj["kz"] = ks_traj["kz"]./(maximum([abs(minimum(ks_traj["kz"])),abs(maximum(ks_traj["kz"]))])*2)
ks_traj = hcat(ks_traj["kx"][:],ks_traj["ky"][:],ks_traj["kz"][:])
ks_traj = permutedims(ks_traj,[2,1])
ks_traj = convert(AbstractMatrix{Float32},ks_traj)
times = repeat(params_pulseq["gen"]["t_vector"][:],numPar)
times = convert(Vector{Float32},times)
ks_traj = Trajectory(ks_traj,1,Int(params_pulseq["gen"]["ro_samples"]*params_pulseq["spi"]["interl"]); 
        times=times,TE=params_pulseq["gen"]["TE"],AQ=params_pulseq["gen"]["ro_time"], 
        numSlices=Int(params_pulseq["gen"]["n_ov"][3]),circular=true)

# Setting repetition label
rawData = SetRepetitionLabel(rawData,params)

# read ME-GRE MRD data and recon
file = ISMRMRDFile(string(params[:path],"/ismrmd/3d/b0_",params[:scan_b0],"_fieldmap.h5"))
acqData_b0 = AcquisitionData(file)
rawData_b0 = RawAcquisitionData(file)
# rawData_b0 = matread(string(params[:path],"/raw/b0_",params[:scan_b0],"_fieldmap.mat"))
# rawData_b0 = rawData_b0["b0"]

# if params_pulseq["gen"]["field_strength"] == 7
#     rawData_b0 = fft(rawData_b0,3)
# end

if params_pulseq["gen"]["field_strength"] == 7
    recon_b0 = ReconCartesianData(acqData_b0,2; interleaved=true)
elseif params_pulseq["gen"]["field_strength"] == 7im
    recon_b0 = ReconCartesianData(acqData_b0,3)
end

# recon_b0 = ReconCartesianData(rawData_b0; dims=2)

# File names
file_name = string(params[:scan_b0],"_",Int(params_pulseq["gen"]["n_ov"][1]),"_",Int(params_pulseq["gen"]["n_ov"][2]),"_",Int(params_pulseq["gen"]["n_ov"][3]))
save_suffix = "" 

####### Option 1) We create the sensitivities with MRIReco
# # Sensitivity maps
# if params[:do_pi_recon]
#     save_suffix = string(save_suffix,"_cs_mrireco")
#     if isfile(string(params[:path],"/acq/cs_",file_name,".jld"))
#         @info ("Sensitivity map exists... Re-calculate? y/n")
#         input = readline()
#     else
#         input = ""
#     end
#     if input == "n"
#         @info ("Loading Sensitivity maps ...")
#         SensitivityMap = load(string(params[:path],"/acq/cs_",file_name,".jld"))
#         SensitivityMap = SensitivityMap["SensitivityMap"]
#     else
#         @info("Calculating Sensitivity Maps...")
#         SensitivityMap = CalculateSensitivityMap(recon_b0,Tuple(Int.(params_pulseq["gen"]["n_ov"])))
#         save(string(params[:path],"/acq/cs_",file_name,".jld"),"SensitivityMap",SensitivityMap)
#     end
# end

### Option 2) We read the sensitivities from Matlab
if params[:do_pi_recon]
    save_suffix = string(save_suffix,"_cs")
    if isfile(string(params[:path],"/acq/cs_",file_name,".mat"))
        @info ("Loading Sensitivity maps ...")
        SensitivityMap = matread(string(params[:path],"/acq/cs_",file_name,".mat"))
        SensitivityMap = SensitivityMap["coil_sens"]
        SensitivityMap = convert(Array{ComplexF32},SensitivityMap)
    else
        error("Sensitivity maps haven't been created... get them and then run script again")
    end
end

# Calculate B0 map
if params[:do_b0_corr]
    save_suffix = string(save_suffix,"_b0")
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
        OffResonanceMap = CalculateOffresonanceMap(recon_b0,SensitivityMap,rawData_b0.params["TE"])
        save(string(params[:path],"/acq/b0_",file_name,".jld"),"OffResonanceMap",OffResonanceMap)
    end
end

# Reconstruction parameters
params_recon = Dict{Symbol, Any}()
params_recon = merge(defaultRecoParams(), params_recon)
if params[:do_pi_recon]
    params_recon[:reco] = "multiCoil"
    params_recon[:senseMaps] = SensitivityMap
end
if params[:do_b0_corr]
    params_recon[:correctionMap] = OffResonanceMap
end
params_pulseq["gen"]["n_ov"] = Int.(params_pulseq["gen"]["n_ov"])
params_recon[:reconSize] = (params_pulseq["gen"]["n_ov"][1],params_pulseq["gen"]["n_ov"][2],params_pulseq["gen"]["n_ov"][3])
params_recon[:regularization] = "L2";
params_recon[:λ] = 1.e-2;  # 1.e-2
params_recon[:iterations] = 20; # (10)
params_recon[:solver] = "cgnr"; # cgnr (L2), fista (L1), admm(L1)
params_recon[:method] = "nfft"; # nfft, leastsquare

# DORK correction
if params[:do_rDORK_corr]
    save_suffix = string(save_suffix,"_rDORK")
    PartRef = 2         # DORK reference repetition
    nav_ref = FormatRawData(rawData,params;single_rep=true,rep_format=PartRef)
    if params[:kz_enc] == 0  # Linear
        nav_ref = nav_ref[6:30,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
        nav_ref = mean(nav_ref, dims=2)
    elseif params[:kz_enc] == 1 # Center-out
        nav_ref = nav_ref[6:30,1,:,1]
        nav_ref = mean(nav_ref, dims=2)
    end
end

# Reconstruct by repetition....
kdata = Array{Array{ComplexF32,2},3}(undef,1,1,1)
if params[:rep_recon] == 0:0 || params[:rep_recon] == 0 
    params[:rep_recon] = 1:params[:numRep]
end
@time begin
    for i_rep=params[:rep_recon]
        tmp = FormatRawData(rawData,params;single_rep=true,rep_format=i_rep)

        # DORK correction
        if params[:do_rDORK_corr]
            tmp = CorrectPartitionDORK(tmp,nav_ref,params)
        end

        kdata[1] = dropdims(reshape(tmp,(:,numCha,1)),dims=3)

        # Create AcqData object
        acqData = AcquisitionData(ks_traj,kdata; encodingSize=(Int(params_pulseq["gen"]["ro_samples"]*params_pulseq["spi"]["interl"]),1,Int(params_pulseq["gen"]["n_ov"][3])), fov=Tuple(params_pulseq["gen"]["fov"]))

        @info("Stop... Before recon...")
        @infiltrate

        # Reconstruction
        @info(string("Reconstructing Repetition # ",i_rep))
        @time Ireco = reconstruction(acqData, params_recon)

        # Lets save the nifti rep by rep
        if params[:do_pi_recon]
            Ireco_mag = NIVolume(abs.(Ireco[:,:,:,1,1,:]))
        else
            Ireco_mag = sqrt.(sum((abs.(Ireco).^2),dims=5))
            Ireco_mag = NIVolume(abs.(Ireco_mag[:,:,:,1,:,:]))
        end

        niwrite(string(params[:path],"/recon/3d/",params[:scan],"_rep_",i_rep,save_suffix,".nii"),Ireco_mag);

        @info("Stop...")
        @infiltrate

        # Trying to clear some memory
        ccall(:malloc_trim, Cvoid, (Cint,), 0) 
        GC.gc()

        ## Save as MAT
        # matwrite(string(params[:path],"/tmp/recon.mat"),Dict("recon" => abs.(Ireco)))
    end
    
    # Merge all repetitions into one file
    TimeSeries = MergeReconstrucion(params[:path],params[:scan],params[:rep_recon],params[:do_pi_recon],params[:do_b0_corr])
    file = string(params[:path],"/recon/",params[:scan],save_suffix,".nii")
    if isfile(file)
        @info ("TimeSeries exists... Re-write? y/n")
        input = readline()
        if input == "y"
            niwrite(file,NIVolume(TimeSeries)) 
        end
    else
        niwrite(file,NIVolume(TimeSeries))
    end

end










# #####################

# @time begin
#     for i=1:length(scans)

#         params[:scan] = scans[i];

#         params[:path] = string(path_tmp,params[:directory]);
#         twix_params = matread(string(params[:path],"acq/",params[:scan],"_twix_params.mat")); twix_params = twix_params["twix_params"];
#         twix_params_b0 = matread(string(params[:path],"acq/",params[:scan_b0],"_twix_params_b0.mat")); twix_params_b0 = twix_params_b0["twix_params_b0"];

#         if params[:traj_type] == "sk"
#             params_pulseq = matread(string(params[:path],"acq/",params[:scan],"_params_sk.mat")); params_pulseq = params_pulseq["params"];
#         else
#             params_pulseq = matread(string(params[:path],"acq/",params[:scan],"_params.mat")); params_pulseq = params_pulseq["params"];
#         end

#         # Some calculations
#         params[:mtx_s] = params_pulseq["gen"]["n_ov"];
#         params[:nx] = Int(params[:mtx_s][1]);
#         params[:ny] = Int(params[:mtx_s][2]);
#         params[:sl] = Int(params[:mtx_s][3].*params_pulseq["gen"]["kz"]) # number of slices
#         params[:sl_no_ov] = Int(params_pulseq["gen"]["n"][3].*params_pulseq["gen"]["kz"])
#         params[:nCoils] = Int(twix_params["ch"]);                   # number of receiver coils
#         params[:kz] = Int(params_pulseq["gen"]["kz"]);
#         params[:repetitions]  = Int(twix_params["repetitions"])    # Number of repetitions
#         params[:dt] = twix_params["dwell"]                          # acquisition dwell time [s]
#         params[:TE] = params_pulseq["gen"]["TE"]; 
#         params[:times] = params_pulseq["gen"]["t_vector"]		    # Times vector for B0 correction
#         params[:seq] = params_pulseq["gen"]["seq"]
#         params[:b0_te] = twix_params_b0["TE"][1:5].*1e-3
#         params[:fov] = params_pulseq["gen"]["fov"]
#         params[:ro_samples] = params_pulseq["gen"]["ro_samples"]
#         params[:field_strength] = params_pulseq["gen"]["field_strength"]
#         params[:echos] = params_pulseq["gen"]["echos"]
#         params[:set] = params_pulseq["gen"]["adc_segments"]
#         # params[:b0_te] = twix_params_b0["TE"][1:9].*1e-3
        
#         # Temp: just for data 0914, to make sure sl are ok
#         if params[:scan_b0] == "s03" && params[:directory] == "data/09142023_sb_9T/"
#             params[:sl] = 142
#         end

#         if params[:seq] == 2
#             params[:contrasts] = ["abc"];
#         end
#         if params[:plt]
#             using Plots, ImageView
#         end
#         if params[:is2d]
#             params[:id] = "2d";
#         else 
#             params[:id] = "3d";
#         end

#         # params[:scan] = string(params[:scan],"_e",Int(2))
#         if params[:echos] > 1
#             scan_orig = params[:scan]
#             for i_ech=1:params[:echos]
#                 params[:scan] = string(scan_orig,"_e",Int(i_ech))
#                 @time fn_sv_recon(params)
#             end
#         else
#             @time fn_sv_recon(params)
#         end

#         ccall(:malloc_trim, Cvoid, (Cint,), 0) 
#         GC.gc()

#     end
# end

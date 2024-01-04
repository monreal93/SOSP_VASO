using Pkg

Pkg.activate("/usr/share/sosp_vaso/recon/")
# cd("/usr/share/5T3/Alejandro/sosp_vaso/")
cd("/usr/share/5T4/Alejandro/sosp_vaso/")
Pkg.activate("./recon/")

using MRIReco, MRIFiles

using MAT, NIfTI #, MriResearchTools
using MRICoilSensitivities
using Revise
using Infiltrator
using MRIFieldmaps: b0map, b0model, b0init, b0scale
using StatsBase: mean
using Unitful: s, ustrip
using MIRTjim: jim, prompt; jim(:prompt, true)
using ImageTransformations
using ImageFiltering: imfilter, Kernel
using NFFT
using FLoops
# using MriResearchTools: makehomogeneous

using FFTW
# FFTW.set_provider!("mkl")
FFTW.set_provider!("fftw")

# This packages are used for the motion correction
# using ImageCore, CoordinateTransformations, Rotations

# NFFT._use_threads[] = false

include("./functions/fn_sv_recon.jl")
include("./functions/fn_save_nii.jl")
include("./functions/fn_prepare_mrd.jl")
include("./functions/fn_calculateSphericalHarmonics.jl")
include("../recon/functions/fn_motionCorrection.jl")
# include("./functions/fn_ismrmd.jl")

params = Dict{Symbol, Any}()
params[:plt] = false
params[:do_pi_recon] = true             # Perform PI reconstruction or direct recon
params[:do_b0_corr] = true
params[:do_t2s_corr] = false
params[:b0_type] = "fessler"               # B0 from: "fessler", "romeo" , "gilad", "skope"               
params[:is2d] = false
# params[:multiRepetitions] = true;         # Reconstruct multiple repetitions, if false = 2nd rep will be reconstructed
params[:rep_recon] = 1                # Range of rep to recon
params[:contrasts] = ["b"]                  # Contrasts to recon v,b,abc or all
params[:traj_type] = "nom"                 # Trajectory type nominal="nom",skope="sk",poet = "poet", corrected = "nom_corr"
params[:save_ph] = 0                       # Save phase of recon as nifti
params[:fmri] = 1                          # 1 for fMRI data will use separate cs and b0 maps per scan
params[:pdork] = ""                         # partition DORK "_pDORK" or ""
params[:rdork] = ""                   # repetition DORK "_rDORK" or ""
params[:idork] = ""                   # interleaves DORK "_iDORK" or ""n_ov
params[:drift] = ""                   # DRIFT correction to B0 map? "_drift" or ""            
params[:mcorr] = ""           # Motion correction with navigators "_mCorr"

# Some parameters
scans = ["sb_01"]            # For now: if multipe echos, include _e1.. _e2..
params[:fieldmap] = "a01"           # Name of the ME-GRE to use for CS and B0map
# scans = ["abc_01","abc_02","abc_03","abc_04","abc_05","abc_06","abc_07","abc_08","abc_09","abc_10","abc_11","abc_12"]; #,"sv_02","sv_04","sv_06","sv_08","sv_09","sv_10"];
params[:directory] = "data/12062023_sk_abc/"        # directory where the data is stored

# Find out if script is running in laptop/dabeast/docker

path_tmp = pwd();
if contains(path_tmp,"amonreal")        # We are in my local laptop
    using Revise, Infiltrator
    path_tmp = SubString(path_tmp,1,38);
elseif contains(path_tmp,"/mnt/")       # We are in daBeast
    params[:plt] = false;
    path_tmp = SubString(path_tmp,1,19);
elseif contains(path_tmp,"/usr/share/")
     # We are in a Docker container
    params[:plt] = false;
    # path_tmp = SubString(path_tmp,1,21);
    path_tmp = string(path_tmp,"/")
end

@time begin
    for i=1:length(scans)

        params[:scan] = scans[i];

        params[:path] = string(path_tmp,params[:directory]);
        twix_params = matread(string(params[:path],"acq/",params[:scan],"_twix_params.mat")); twix_params = twix_params["twix_params"];
        twix_params_b0 = matread(string(params[:path],"acq/",params[:fieldmap],"_twix_params_b0.mat")); twix_params_b0 = twix_params_b0["twix_params_b0"];

        if params[:traj_type] == "sk"
            pulseq_params = matread(string(params[:path],"acq/",params[:scan],"_params_sk.mat")); pulseq_params = pulseq_params["params"];
        else
            pulseq_params = matread(string(params[:path],"acq/",params[:scan],"_params.mat")); pulseq_params = pulseq_params["params"];
        end

        # Some calculations
        params[:mtx_s] = pulseq_params["gen"]["n_ov"];
        params[:nx] = Int(params[:mtx_s][1]);
        params[:ny] = Int(params[:mtx_s][2]);
        params[:sl] = Int(params[:mtx_s][3].*pulseq_params["gen"]["kz"]) # number of slices
        params[:sl_no_ov] = Int(pulseq_params["gen"]["n"][3].*pulseq_params["gen"]["kz"])
        params[:nCoils] = Int(twix_params["ch"]);                   # number of receiver coils
        params[:kz] = Int(pulseq_params["gen"]["kz"]);
        params[:repetitions]  = Int(twix_params["repetitions"])    # Number of repetitions
        params[:dt] = twix_params["dwell"]                          # acquisition dwell time [s]
        params[:TE] = pulseq_params["gen"]["TE"]; 
        params[:times] = pulseq_params["gen"]["t_vector"]		    # Times vector for B0 correction
        params[:seq] = pulseq_params["gen"]["seq"]
        params[:b0_te] = twix_params_b0["TE"][1:5].*1e-3
        params[:fov] = pulseq_params["gen"]["fov"]
        params[:ro_samples] = pulseq_params["gen"]["ro_samples"]
        params[:field_strength] = pulseq_params["gen"]["field_strength"]
        params[:echos] = pulseq_params["gen"]["echos"]
        # params[:b0_te] = twix_params_b0["TE"][1:9].*1e-3
        
        # Temp: just for data 0914, to make sure sl are ok
        if params[:fieldmap] == "s03" && params[:directory] == "data/09142023_sb_9T/"
            params[:sl] = 142
        end

        if params[:seq] == 2
            params[:contrasts] = ["abc"];
        end
        if params[:plt]
            using Plots, ImageView
        end
        if params[:is2d]
            params[:id] = "2d";
        else 
            params[:id] = "3d";
        end

        # params[:scan] = string(params[:scan],"_e",Int(2))
        if params[:echos] > 1
            scan_orig = params[:scan]
            for i_ech=1:params[:echos]
                params[:scan] = string(scan_orig,"_e",Int(i_ech))
                @time fn_sv_recon(params)
            end
        else
            @time fn_sv_recon(params)
        end

        ccall(:malloc_trim, Cvoid, (Cint,), 0) 
        GC.gc()

    end
end

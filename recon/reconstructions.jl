
using MRIReco, MAT, NIfTI, MriResearchTools
using Revise
using Infiltrator
using MRIFieldmaps: b0map, b0model, b0init, b0scale
using StatsBase: mean
using Unitful: s, ustrip
using MIRTjim: jim, prompt; jim(:prompt, true)
# using MIRT: ir_mri_sensemap_sim

using MRIFiles

include("./functions/fn_sv_recon.jl")
include("./functions/fn_save_nii.jl")
# include("./functions/fn_ismrmd.jl")

params = Dict{Symbol, Any}()
params[:plt] = false;
params[:do_b0_corr] = false;
params[:b0_type] = "fessler";                 # B0 from: "fessler", "romeo" , "gilad", "skope"               
params[:is2d] = false;
params[:multiRepetitions] = false;          # Reconstruct multiple repetitions, if false = 2nd rep will be reconstructed
params[:rep_recon] = 1;                    # Rep to recon if multiRepetitions=false
params[:contrasts] = ["v"];                 # Contrasts to recon v,b or both 
params[:traj_type] = "nom";                 # Trajectory type nominal="nom",skope="sk"
params[:save_ph] = 0;                       # Save phase of recon as nifti

# Some parameters
scans = ["sv_01"]; #,"sv_02","sv_03","sv_04","sv_05"];
params[:directory] = "data/10292022_sv_2/"        # directory where the data is stored

# Find out if script is running in laptop/dabeast/docker
path_tmp = pwd();
if contains(path_tmp,"amonreal")        # We are in my local laptop
    using Revise, Infiltrator
    path_tmp = SubString(path_tmp,1,38);
elseif contains(path_tmp,"/mnt/")       # We are in daBeast
    params[:plt] = false;
    path_tmp = SubString(path_tmp,1,19);
elseif contains(path_tmp,"/usr/share/") # We are in a Docker container
    params[:plt] = false;
    path_tmp = SubString(path_tmp,1,21);
end

@time begin
    for i=1:length(scans)

        params[:scan] = scans[i];

        params[:path] = string(path_tmp,params[:directory]);
        twix_params = matread(string(params[:path],"acq/",params[:scan],"_twix_params.mat")); twix_params = twix_params["twix_params"];
        pulseq_params = matread(string(params[:path],"acq/",params[:scan],"_params.mat")); pulseq_params = pulseq_params["params"];

        # Some calculations
        mtx_s = pulseq_params["gen"]["n"];  
        params[:nx] = Int(mtx_s[1]);
        params[:ny] = Int(mtx_s[2]);
        params[:sl] = Int(mtx_s[3]); #*pulseq_params["gen"]["kz"]); # number of slices
        params[:nCoils] = Int(twix_params["ch"]);                   # number of receiver coils
        params[:kz] = Int(pulseq_params["gen"]["kz"]);
        params[:repetitions]  = Int(twix_params["repetitions"]);    # Number of repetitions
        params[:dt] = twix_params["dwell"]                          # acquisition dwell time [s]
        params[:TE] = pulseq_params["gen"]["TE"]; 
        params[:times] = pulseq_params["gen"]["t_vector"];		    # Times vector for B0 correction
        params[:seq] = pulseq_params["gen"]["seq"];
        params[:b0_te] = pulseq_params["twix_params_b0"]["TE"][1:9].*1e-3

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

        @time fn_sv_recon(params)

        GC.gc()

    end
end

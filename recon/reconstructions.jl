
using MRIReco, MAT, NIfTI, MriResearchTools
using Revise
using Infiltrator

include("./functions/fn_sv_recon.jl")
include("./functions/fn_save_nii.jl")

params = Dict{Symbol, Any}()

params[:plt] = false;
params[:do_b0_corr] = false;
params[:b0_type] = "romeo";                     # B0 map from: "romeo" , "gilad", "skope"               
params[:is2d] = true;
params[:multiRepetitions] = false;              # Reconstruct multiple repetitions, if false = 2nd rep will be reconstructed

# AMM: Temp: For now I set this param for 2d Recon manually, I should take it from create_ismrmd_cs_b0_v1.m 
params[:sl_reco] = 5;

# Some parameters
params[:scan] = "sv_02";                       # sv_#_b (scan # Bold) or sv_#_v (scan # Vaso)
params[:directory] = "data/sv_05302022/"        # directory where the data is stored

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

params[:path] = string(path_tmp,params[:directory]);
twix_params = matread(string(params[:path],"acq/",params[:scan],"_twix_params.mat")); twix_params = twix_params["twix_params_sv"];
pulseq_params = matread(string(params[:path],"acq/",params[:scan][1:5],"_params.mat")); pulseq_params = pulseq_params["params"];

mtx_s = twix_params["mtx_s"];  
params[:nx] = Int(mtx_s[1]);                                 # sv1 = 112,112 sv9, 218,216
params[:ny] = Int(mtx_s[2]);
params[:sl] = Int(mtx_s[3]);                                 # number of slices
params[:nCoils] = Int(twix_params["ch"]);                    # number of receiver coils
params[:repetitions]  = Int(twix_params["repetitions"]);     # Number of repetitions
# params[:dt] = twix_params["dwell"]                      # acquisition dwell time [s]
params[:dt] = 2e-6;                      		 # acquisition dwell time [s]
params[:TE] = pulseq_params["gen"]["TE"]; 
# params[:sl_reco] = Int(twix_params["sl_to_recon"]);                          # slice to reconstruct if is2d = true
params[:times] = pulseq_params["gen"]["t_vector"];		# Times vector for B0 correction

if params[:plt]
    using Plots, ImageView
end

if params[:is2d]
    params[:id] = "2d";
else 
    params[:id] = "3d";
end

fn_sv_recon(params)

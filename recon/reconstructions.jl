
using MRIReco, MAT, NIfTI, MriResearchTools
using Revise
using Infiltrator

include("./functions/fn_sv_recon.jl")
include("./functions/fn_save_nii.jl")
# include("./functions/fn_ismrmd.jl")

params = Dict{Symbol, Any}()

# @info string("Reconstructing scan ", ARGS[1], ", contrast ", ARGS[2])

params[:plt] = false;
params[:do_b0_corr] = true;
params[:b0_type] = "romeo";                     # B0 map from: "romeo" , "gilad", "skope"               
params[:is2d] = false;
params[:multiRepetitions] = false;              # Reconstruct multiple repetitions, if false = 2nd rep will be reconstructed
params[:contrasts] = ["v"];                # Contrasts to recon v,b or both 
# contrasts = ["b","v"];
# params[:contrasts] = ARGS[2];

# AMM: Temp: For now I set this param for 2d Recon manually, I should take it from create_ismrmd_cs_b0_v1.m 
params[:sl_reco] = 7;

# Some parameters
# params[:scan] = "sv_01";                       # sv_#_b (scan # Bold) or sv_#_v (scan # Vaso)
# scans = ["sv_01","sv_02","sv_03","sv_04","sv_05","sv_06"];
scans = ["sv_04"];#,"sv_02","sv_03","sv_04","sv_05"];
# params[:scan] = ARGS[1];
params[:directory] = "data/08032022/"        # directory where the data is stored

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
        twix_params = matread(string(params[:path],"acq/",params[:scan],"_twix_params.mat")); 
        twix_params = twix_params["twix_params"];
        pulseq_params = matread(string(params[:path],"acq/",params[:scan],"_params.mat")); pulseq_params = pulseq_params["params"];

        # Some calculations
        mtx_s = pulseq_params["gen"]["n"];  
        params[:nx] = Int(mtx_s[1]);                                 # sv1 = 112,112 sv9, 218,216
        params[:ny] = Int(mtx_s[2]);
        params[:sl] = Int(mtx_s[3]); #*pulseq_params["gen"]["kz"]);   # number of slices
        params[:nCoils] = Int(twix_params["ch"]);                    # number of receiver coils
        params[:kz] = Int(pulseq_params["gen"]["kz"]);
        params[:repetitions]  = Int(twix_params["repetitions"]);     # Number of repetitions
        # params[:dt] = twix_params["dwell"]                      # acquisition dwell time [s]
        params[:dt] = 2e-6;                      		 # acquisition dwell time [s]
        params[:TE] = pulseq_params["gen"]["TE"]; 
        # params[:sl_reco] = Int(twix_params["sl_to_recon"]);                          # slice to reconstruct if is2d = true
        params[:times] = pulseq_params["gen"]["t_vector"];		# Times vector for B0 correction
        params[:seq] = pulseq_params["gen"]["seq"];
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

        # Load ISMRMD file and reshape it
        # fn_ismrmd(params)

        # Reconstruction
        # params[:contrasts] = contrasts[i];
        
        @time fn_sv_recon(params)

        GC.gc()

    end
end
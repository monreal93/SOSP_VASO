using Revise
using MRIReco, MAT, NIfTI, MriResearchTools
using Infiltrator
# using Plots, ImageView

multicoil = true;
sl = 10;
is2d = true;

# Load volume for simulation
phn = matread("/usr/share/sosp_vaso/data/sv_05302022/tmp/gre_250_250_30_9_32.mat");
phn = phn["gre"];
phn = convert(Array{ComplexF64,5},phn);

# Load trajectory
# traj = matread("/home/amonreal/Documents/PhD/PhD_2022/data/tmp/sv_01_ks_traj_nom_norm.mat");
traj = matread("/usr/share/sosp_vaso/data/sv_05302022/acq/sv_01_ks_traj_nom.mat");
traj = traj["ks_traj"];
traj["kx"] = traj["kx"]./maximum(traj["kx"])./2;
traj["ky"] = traj["ky"]./maximum(traj["ky"])./2;
traj["kz"] = traj["kz"]./maximum(traj["kz"])./2;

# Load Coil Sensitivities
cs = matread("/usr/share/sosp_vaso/data/sv_05302022/acq/cs_250_250_30.mat");
cs = cs["coil_sens"];

# Load B0 map
b0 = niread("/usr/share/sosp_vaso/data/sv_05302022/acq/b0_250_250_30_gilad.nii");
b0 = b0.raw;

# Load parameters
params_pulseq = matread("/usr/share/sosp_vaso/data/sv_05302022/acq/sv_01_params.mat");
params_pulseq = params_pulseq["params"];

# Subseting if 2D
if is2d
    phn = phn[:,:,sl,1,:];
    phn = sqrt.(sum((phn.^2),dims=3));
    phn = dropdims(phn,dims=3);
    # phn = reshape(phn,(size(phn)...,1));
    # phn = permutedims(phn, (1,2,4,3));
    cs = cs[:,:,sl,:];
    cs = reshape(cs,(size(cs)...,1));
    cs = permutedims(cs, (1,2,4,3));
    b0 = b0[:,:,sl,:];
    # b0 = reshape(b0,(size(b0)...,1));
    # b0 = permutedims(b0, (1,2,4,3));
    b0 = 1im*b0.*2*pi.*4;
    traj["kx"] = traj["kx"][:,sl];
    traj["ky"] = traj["ky"][:,sl];
    traj["kz"] = traj["kz"][:,sl];
end

if is2d
    ks = [vec(traj["kx"]) vec(traj["ky"])];
else
    ks = [vec(traj["kx"]) vec(traj["ky"]) vec(traj["kz"])];
end
ks = permutedims(ks,(2,1));
ks = Trajectory(ks,Int64(params_pulseq["spi"]["interl"]),Int64(params_pulseq["spi"]["ro_samples"]));

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Spiral"
params[:AQ] = 30.0e-3;
# params[:numProfiles] = 1
# params[:numSamplingPerProfile] = N*N
# params[:windings] = 128
# params[:AQ] = 3.0e-2

# do simulation
if multicoil
    acqData = simulation(ks, phn, correctionMap = b0[:,:,1]; senseMaps=cs, params)
else
    acqData = simulation(ks, phn, correctionMap = []; senseMaps=[], params)
end

# reco parameters
params = Dict{Symbol, Any}()

if is2d
    params[:reconSize] = (Int(params_pulseq["gen"]["n"][1]),Int(params_pulseq["gen"]["n"][2]))
else
    params[:reconSize] = (params_pulseq["gen"]["n"][1],params_pulseq["gen"]["n"][2],params_pulseq["gen"]["n"][3])
end
params[:regularization] = "L2"
params[:λ] = 1.e-3
params[:iterations] = 40
params[:solver] = "cgnr"

if multicoil
    params[:reco] = "multiCoil"
    params[:senseMaps] = cs;
    #params[:correctionMap] = b0;
end

# do reconstruction
Ireco = reconstruction(acqData, params)

Ireco_nii = NIVolume(abs.(Ireco[:,:,:,1,1]));

niwrite("/usr/share/sosp_vaso/sim/recon.nii",Ireco_nii)
Pkg.activate("/usr/share/sosp_vaso/sim/")

using Revise, Pkg
using MRIReco, MRISimulation, MAT, NIfTI, ImagePhantoms
using MIRTjim: jim, prompt; jim(:prompt, true)
using Infiltrator
using ImageTransformations
# using Plots, ImageView

include("../recon/functions/fn_addCorrelatedNoise.jl")
include("../recon/functions/fn_calculateGfactor.jl")
include("../recon/functions/fn_save_nii.jl")

phn_sim = 1             # 1=brain, 0=point (psf)
cs_sim = true
cs_recon = true
b0_sim = false
b0_recon = false
add_noise = false
gfactor = false
sl = 10
is2d = false
# Folder and name of sequence to simulate
folder_sim = "zoia"
scan_sim = ["sv_09"]
# scan_sim = ["sv_01","sv_02","sv_03","sv_04","sv_05","sv_06","sv_07","sv_08","sv_09","sv_10","sv_11"]

# Folder and name of sensitivity maps and b0 map to use for simulation
folder = "03032023_sv"
scan = "sv_01"
path = string("/usr/share/sosp_vaso/data/",folder)
path_sim = string("/usr/share/sosp_vaso/data/",folder_sim)

##### Load parameters
params = matread(string(path,"/acq/",scan,"_params.mat"))
params = params["params"]
scan_suffix = string(scan,"_",Int(params["gen"]["n"][1]),"_",Int(params["gen"]["n"][2]),"_",Int(params["gen"]["n"][3]))

@infiltrate
###### Load volume for simulation
if phn_sim == 1
    phn = matread(string(path,"/acq/fm_",scan_suffix,".mat"))
    phn = phn["fieldmap"]
    phn = dropdims(sqrt.(sum(abs.(phn).^2; dims=5));dims=5)
    phn = phn[:,:,:,1]
    phn = convert(Array{ComplexF64,3},phn)
    phn = imresize(phn,(Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3])))
    phn_abs = abs.(phn)
    phn_abs = (phn_abs.- minimum(last,phn_abs))./(maximum(last,phn_abs)-minimum(last,phn_abs))
    phn_nii = NIVolume(phn_abs)
    niwrite(string("/usr/share/sosp_vaso/data/",folder_sim,"/sim/reference_phn.nii"),phn_nii)
else
    phn = zeros(mtx_s)
    phn[Int(params["gen"]["n"][1]/2),Int(params["gen"]["n"][2]/2),Int(params["gen"]["n"][3]/2)] = 1;
    phn = Array{ComplexF64}(phn)
    phn_nii = NIVolume(phn)
    niwrite(string("/usr/share/sosp_vaso/data/",folder_sim,"/sim/reference_phn_point.nii"),phn_nii)
end

###### Load Coil Sensitivities
cs = matread(string(path,"/acq/cs_",scan_suffix,".mat"))
cs = cs["coil_sens"]
cs = convert(Array{ComplexF64,4},cs)
cs = imresize(cs,(Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3])))
# cs = reverse(cs,dims = 2)

##### Load B0 map
if b0_sim
    b0 = matread(string(path,"/acq/b0_",scan_suffix,".mat"))
    b0 = b0["b0"]
    b0 = convert(Array{ComplexF64,3},b0)
end

# Subseting if 2D
if is2d
    phn = phn[:,:,sl,1,:]
    phn = sqrt.(sum((phn.^2),dims=3))
    phn = dropdims(phn,dims=3)
    # phn = reshape(phn,(size(phn)...,1));
    # phn = permutedims(phn, (1,2,4,3));
    cs = cs[:,:,sl,:]
    cs = reshape(cs,(size(cs)...,1))
    cs = permutedims(cs, (1,2,4,3))
    b0 = b0[:,:,sl,:]
    # b0 = reshape(b0,(size(b0)...,1));
    # b0 = permutedims(b0, (1,2,4,3));
    b0 = 1im*b0.*2*pi.*4
    traj["kx"] = traj["kx"][:,sl]
    traj["ky"] = traj["ky"][:,sl]
    traj["kz"] = traj["kz"][:,sl]
end

# if is2d
#     ks = [vec(traj["kx"]) vec(traj["ky"])]
# else
#     ks = [vec(traj["kx"]) vec(traj["ky"]) vec(traj["kz"])]
# end

for i = 1:Int(length(scan_sim))

    ###### Load simulation parameters
    params_pulseq = matread(string(path_sim,"/acq/",scan_sim[i],"_params.mat"))
    params_pulseq = params_pulseq["params"]

    ###### Load trajectory
    traj = matread(string(path_sim,"/acq/",scan_sim[i],"_ks_traj_nom.mat"))
    traj = traj["ks_traj"]

    # Normalizing the traj to -0.5 to 0.5
    traj = [traj["kx"][:] traj["ky"][:] traj["kz"][:]]
    traj = permutedims(traj,[2,1])

    kx_fov = maximum([abs.(minimum(traj[1,:])) abs.(maximum(traj[1,:]))])*2
    ky_fov = maximum([abs.(minimum(traj[2,:])) abs.(maximum(traj[2,:]))])*2
    kz_fov = maximum([abs.(minimum(traj[3,:])) abs.(maximum(traj[3,:]))])*2

    traj[1,:] = traj[1,:]./kx_fov 
    traj[2,:] = traj[2,:]./ky_fov 
    traj[3,:] = traj[3,:]./kz_fov

    # Creating/loading time vector
    times = params_pulseq["gen"]["t_vector"][:]		    # Times vector for B0 correction
    # n_samples = Int(size(traj)[2]/params["gen"]["n"][3]/params["spi"]["interl"])
    # tAQ = (n_samples-1) * params["gen"]["adc_dwell"]
    # times = params["gen"]["te"] .+ collect(0:params["gen"]["adc_dwell"]:tAQ);

    # if its 3D, repeat the times vector for each slice
    if !is2d
        times = repeat(times,Int(params_pulseq["gen"]["n"][3]/params_pulseq["gen"]["kz"]))
        times = vec(times)
    end

    ks = permutedims(traj,(2,1))
    ks = Trajectory(traj,Int64(params_pulseq["spi"]["interl"]),Int64(params_pulseq["gen"]["ro_samples"]);times = times)

    ###### Simulation
    params_sim = Dict{Symbol, Any}()
    params_sim[:simulation] = "fast"
    params_sim[:trajName] = "Spiral"
    params_sim[:AQ] = params["gen"]["acqTR"]
    params_sim[:times] = times

    # Do simulation
    @info ("Simulating k-space data ...")
    senseMaps = []; correctionMap = [];
    if cs_sim
        # acqData = simulation(ks, phn, correctionMap = []; opName="fast", senseMaps=cs, params_sim);
        # acqData = simulation(ks, phn, correctionMap = b0[:,:,1]; senseMaps=cs, params)
        senseMaps = cs
    end
    if b0_sim
        correctionMap = b0
        acqData = simulation(ks, phn, correctionMap; senseMaps=senseMaps, params_sim)
    else
        acqData = simulation(ks, phn, correctionMap=[]; senseMaps=senseMaps, params_sim)
    end

    @infiltrate

    if params_pulseq["gen"]["ro_type"] == "s"
        noise_snr = Float64(5*(params_pulseq["spi"]["interl"]))
        if params_pulseq["spi"]["interl"] > 1
            noise_snr = noise_snr*2
        end
    elseif params_pulseq["gen"]["ro_type"] == "c"
        noise_snr = Float64(2)
    end

    ##### Add Correlated noise
    if add_noise
        # load covariance Matrix
        cov = matread(string("/usr/share/sosp_vaso/data/tmp/covariance.mat")); cov = cov["C"]
        # acqData = addNoise(acqData,noise_snr)
        acqData.kdata[1] = addCorrelatedNoise(acqData.kdata[1],noise_snr,cov,Float64(1))
    end

    ###### Reconstruction
    params_reco = Dict{Symbol, Any}()

    if is2d
        params_reco[:reconSize] = Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2])
    else
        params_reco[:reconSize] = Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3])
    end
    params_reco[:regularization] = "L2"
    params_reco[:λ] = 1.e-3
    params_reco[:iterations] = 40
    params_reco[:solver] = "cgnr"
    params_reco[:method] = "nfft"

    if cs_recon
        params_reco[:reco] = "multiCoil"
        params_reco[:senseMaps] = cs
    end

    if b0_recon
        params_reco[:correctionMap] = b0
    end

    # G-factor map
    if gfactor
        @info ("Calculating G-factor ...")
        g_factor = calculateGfactor(acqData,6,params_reco)
    end

    # Do reconstruction
    @info ("Reconstruction ...")
    Ireco = reconstruction(acqData, params_reco)
    @infiltrate

    ###### Quality measurments
    # Normalize datasets
    if cs_recon
        Ireco = abs.(Ireco[:,:,:,1,1])
    else
        Ireco = abs.(Ireco[:,:,:,1,:,1])
    end
    Ireco = (Ireco.- minimum(last,Ireco))./(maximum(last,Ireco)-minimum(last,Ireco))
    recon_error = Ireco.-phn_abs

    suffix = "_recon"
    if phn_sim == 0
        suffix = string(suffix,"_point")
    end
    if cs_sim
        suffix = string(suffix,"_cssim")
    end
    if cs_recon
        suffix = string(suffix,"_csrecon")
    end
    if b0_sim
        suffix = string(suffix,"_b0sim")
    end
    if b0_recon
        suffix = string(suffix,"_b0recon")
    end

    # G-factor map
    if gfactor
        g_factor = calculateGfactor(acqData,6,params_reco)
    end

    if cs_sim == true && cs_recon ==false
        Ireco = sqrt.(sum((Ireco.^2),dims=4))
        Ireco = Ireco[:,:,:,1]
    end

    @infiltrate

    # save recon
    Ireco_nii = NIVolume(abs.(Ireco))
    path_save_recon = string("/usr/share/sosp_vaso/data/",folder_sim,"/sim/",scan_sim[i],"_sim",suffix,".nii")
    niwrite(path_save_recon,Ireco_nii)
    
    # save error
    recon_error = NIVolume(recon_error);
    path_save_error = string("/usr/share/sosp_vaso/data/",folder_sim,"/sim/",scan_sim[i],"_error",suffix,".nii")
    niwrite(path_save_error,recon_error)
end
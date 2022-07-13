
using MRIReco, MAT, NIfTI

N = 256
I = shepp_logan(N)
I = circularShutterFreq!(I,1)
cmap = 1im*quadraticFieldmap(N,N,125*2pi)

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Spiral"
params[:numProfiles] = 1
params[:numSamplingPerProfile] = N*N
params[:windings] = 128
params[:AQ] = 3.0e-2
params[:correctionMap] = cmap[:,:,1]

# do simulation
acqData = simulation(I, params)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N,N)
params[:correctionMap] = cmap
params[:alpha] = 1.75
params[:m] = 4.0
params[:K] = 28

# Temp: AMM: Messing up with the time vector
acqData.traj[1].times = acqData.traj[1].times.*1e3;

# do reconstruction
Ireco = reconstruction(acqData, params)

Ireco_nii = NIVolume(abs.(Ireco[:,:,:,1,1]));

niwrite("/usr/share/sosp_vaso/sim/mrireco_recon.nii",Ireco_nii)
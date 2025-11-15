include("fn_addCorrelatedNoise.jl")

using Infiltrator

"""
  Calculates g-factor maps using the Pseudo Multiple Replica method

# Arguments
* `acqData::AcquisitionData`    - AcquisitionData
* 'params::Dict'          - reconstruction parameters

"""
function calculateGfactor(acqData::AcquisitionData,acqData_full::AcquisitionData, replicas::Int64, cov::Matrix, params::Dict)

  image_replicas = Array{ComplexF64}(undef,params[:reconSize][1],params[:reconSize][2],params[:reconSize][3],replicas)
  noise_replicas = Array{ComplexF64}(undef,params[:reconSize][1],params[:reconSize][2],params[:reconSize][3],replicas)
  
  noise_snr = Float64(10) # SNR for simulations...
 
  # Reconstruct image data
  # Ireco_no_noise = reconstruction(acqData, params)

  # Reconstruct noise data
  # file_name = "/neurodesktop-storage/5T4/Alejandro/sosp_vaso/data/simulations/tmp/"
  idx = findall.("_",scan_suffix)
  if params[:ro_type] == "s"
    file_name = string(path_sim,"/tmp/spiral_noise_replicas",scan_suffix[idx[1][1]:end],".mat")
  elseif params[:ro_type] == "c"
    file_name = string(path_sim,"/tmp/cartesian_noise_replicas",scan_suffix[idx[1][1]:end],".mat")
  end

  signalAmpl = sum(abs.(acqData.kdata[1])) / length(acqData.kdata[1])

  if !isfile(file_name)
    acqData_zeros = deepcopy(acqData_full)
    acqData_zeros.kdata[1] .= zeros(size(acqData_zeros.kdata[1]))
    acqData_zeros_noise = deepcopy(acqData_zeros)
    
    for i_rep=1:replicas
    # @floop for i_rep=1:replicas
      scale_factor = 1e-3

      acqData_zeros_noise.kdata[1] .= addCorrelatedNoise(acqData_zeros.kdata[1],noise_snr,cov; signalAmpl=signalAmpl, scale_factor=scale_factor)

      Ireco = reconstruction(acqData_zeros_noise, params)
      
      Ireco *= sqrt(params[:rxyz])

      noise_replicas[:,:,:,i_rep] = (Ireco[:,:,:,1,1,1])
      @info (string("Done with noise only replica ",i_rep))

    end
    
    # Temp: Saving noise_replicas 
    matwrite(file_name, Dict("noise_replicas" => noise_replicas))
  else
    noise_replicas = matread(file_name)
    noise_replicas = noise_replicas["noise_replicas"]
    
  end
 
  @info("stop g-factor....")
  @infiltrate

  # Reconstruct image-noise data
  acqData_noise = deepcopy(acqData)
  for i_rep=1:replicas
  # @floop for i_rep=1:replicas
      scale_factor = 1e-3
      acqData_noise.kdata[1] .= addCorrelatedNoise(acqData.kdata[1],noise_snr,cov;  scale_factor=scale_factor)

      Ireco = reconstruction(acqData_noise, params)

      image_replicas[:,:,:,i_rep] = (Ireco[:,:,:,1,1,1].data)
      @info (string("Done with image-noise replica ",i_rep))

  end

  @info("stop g-factor....")
  @infiltrate

  # Get standard deviation of image-noise and noise reconstructions
  noise_std = std(real.(noise_replicas); dims=4)
  image_std = std(real.(image_replicas); dims=4)

  gfactor = image_std./(noise_std.*sqrt(params[:rxyz]))
  gfactor[isnan.(gfactor)] .= 0
  gfactor = gfactor[:,:,:,1]

  return gfactor
end
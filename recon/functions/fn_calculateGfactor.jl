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
  noise_snr = Float64(5)

 
  # Reconstruct image data
  # Ireco_no_noise = reconstruction(acqData, params)

  # Reconstruct noise data
  file_name = "/usr/share/5T4/Alejandro/sosp_vaso/data/simulations/tmp/"
  if params[:scan_suffix][1] == 's'
    file_name = string(params[:path_sim],"/tmp/spiral_noise_replicas",params[:scan_suffix][6:end],".mat")
  elseif params[:scan_suffix][1] == 'c'
    file_name = string(params[:path_sim],"/tmp/cartesian_noise_replicas",params[:scan_suffix][6:end],".mat")
  end
  if !isfile(file_name)
    acqData_zeros = deepcopy(acqData_full)
    acqData_zeros.kdata[1] .= zeros(size(acqData_zeros.kdata[1]))
    acqData_zeros_noise = deepcopy(acqData_zeros)
    for i_rep=1:replicas
      scale_factor = 1e1;
      acqData_zeros_noise.kdata[1] .= addCorrelatedNoise(acqData_zeros.kdata[1],noise_snr,cov,scale_factor)

      Ireco = reconstruction(acqData_zeros_noise, params)
      
      Ireco *= sqrt(params[:rxyz])

      noise_replicas[:,:,:,i_rep] = (Ireco[:,:,:,1,1,1])
      # acqData_zeros_noise.kdata[1] .= acqData_zeros.kdata[1]
      @info (string("Done with noise only replica ",i_rep))

    end
    
    # Temp: Saving noise_replicas 
    matwrite(file_name, Dict("noise_replicas" => noise_replicas))
  else
    noise_replicas = matread(file_name)
    noise_replicas = noise_replicas["noise_replicas"]
  end

  # Reconstruct image-noise data
  acqData_noise = deepcopy(acqData)
  for i_rep=1:replicas
      scale_factor = 1e1;
      acqData_noise.kdata[1] .= addCorrelatedNoise(acqData.kdata[1],noise_snr,cov,scale_factor)

      Ireco = reconstruction(acqData_noise, params)

      image_replicas[:,:,:,i_rep] = (Ireco[:,:,:,1,1,1].data)
      @info (string("Done with image-noise replica ",i_rep))

  end

  @info("stop....")
  @infiltrate

  # Get standard deviation of image-noise and noise reconstructions
  noise_std = std(real.(noise_replicas); dims=4)
  image_std = std(real.(image_replicas); dims=4)

  gfactor = image_std./(noise_std.*sqrt(params[:rxyz]))
  gfactor[isnan.(gfactor)] .= 0
  gfactor = gfactor[:,:,:,1]

  return gfactor
end
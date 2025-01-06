using LinearAlgebra
using Infiltrator

"""
  Adds correlated average white gaussian noise to the signal x

# Arguments
* `x::Vector`     - signal vector
* 'snr::Float64'  - target SNR

"""
function addCorrelatedNoise(x::Matrix, snr::Float64 ,cov::Matrix; signalAmpl = nothing, scale_factor=1, complex= true)

  if isnothing(signalAmpl)
    signalAmpl = sum(abs.(x))/length(x)
  end

  @info("Stop.. addCorrelatedNoise..")
  @infiltrate

  # target noise amplitude
  noiseAmpl = signalAmpl/snr

  if complex
    # noise = noiseAmpl/sqrt(2.)*( randn(size(x))+ 1im*randn(size(x)) )
    noise = rand(size(x)[1],size(x)[2])+ (1im*rand(size(x)[1],size(x)[2]))
  else
    noise = rand(size(x)[1],size(x)[2])
  end
  
  ## option 1) correlating noise, using cov to scale data, for this option scale factor should be 1e-3,
  ## to match the signal amplitude of the simulated data 
  eig = eigen(cov)
  cov_sqrt = eig.vectors * sqrt.(Diagonal(eig.values)) * inv(eig.vectors)
  noise = permutedims(noise,(2,1))
  noise = cov_sqrt*noise
  noise = permutedims(noise,(2,1))
  noise = noise.*scale_factor

  # ## option 2) Just off-diagonal elements.. and scaling using SNR value
  # noise .*= noiseAmpl
  # cov[diagind(cov)] .= 1
  # noise = permutedims(noise,(2,1))
  # noise = cov*noise
  # noise = permutedims(noise,(2,1))
  # noise = noise.*scale_factor

  return x.+noise
end
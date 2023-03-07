using Infiltrator

"""
  Adds correlated average white gaussian noise to the signal x

# Arguments
* `x::Vector`     - signal vector
* 'snr::Float64'  - target SNR

"""
function addCorrelatedNoise(x::Matrix, snr::Float64 ,cov::Matrix, scale_factor::Float64, complex= true)
  signalAmpl = sum(abs.(x))/length(x)

  # target noise amplitude
  noiseAmpl = signalAmpl/snr

  if complex
    # noise = noiseAmpl/sqrt(2.)*( randn(size(x))+ 1im*randn(size(x)) )
    noise = randn(size(x))+ 1im*randn(size(x))
  else
    noise = noiseAmpl*randn(size(x))
  end

  # correlating noise
  noise = permutedims(noise,(2,1))
  noise = sqrt(cov)*noise
  noise = permutedims(noise,(2,1))
  noise = noise.*scale_factor
  
  @infiltrate

  return x+noise
end
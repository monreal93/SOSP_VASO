include("fn_addCorrelatedNoise.jl")

using Infiltrator

"""
  Calculates g-factor maps using the Pseudo Multiple Replica method

# Arguments
* `acqData::AcquisitionData`    - AcquisitionData
* 'params::Dict'          - reconstruction parameters

"""
function calculateGfactor(acqData::AcquisitionData, replicas::Int64, params::Dict)

    image_replicas = Array{Float64}(undef,params[:reconSize][1],params[:reconSize][2],params[:reconSize][3],replicas)
    noise_snr = Float64(5)
    for i_rep=1:replicas
        scale_factor = sqrt(2).*i_rep
        cov = matread(string("/usr/share/sosp_vaso/data/tmp/covariance.mat")); cov = cov["C"]
        acqData.kdata[1] = addCorrelatedNoise(acqData.kdata[1],noise_snr,cov,scale_factor)

        @info ("Calculating G-factor map ...")
        Ireco = reconstruction(acqData, params)
        image_replicas[:,:,:,i_rep] = abs.(Ireco[:,:,:,1,1,1])

        @infiltrate

    end

end
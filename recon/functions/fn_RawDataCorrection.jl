
"""
CorrectFFTShift(RawData::RawAcquisitionData,k0::AbstractVector,params::Dict{Symbol,Any})

Performs FFT shift
"""

function CorrectFFTShift(tmp,fftShift::Tuple,ksTraj::Matrix{Float32},params::Dict{Symbol,Any})
    
    ksTraj = ksTraj[:,1:(params[:numRead]*params[:numInterl])]'

    tmp = tmp.*exp.(1im*2*π.*ksTraj[:,1].*(fftShift[2]))  # kx shift
    tmp = tmp.*exp.(1im*2*π.*ksTraj[:,2].*(fftShift[1]))  # ky shift
    # tmp = tmp.*exp.(1im*2*π.*ksTraj[3,:].*fftShift[3])  # kz shift
    
    return tmp
end


"""
CorrectRepetitionDORK(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function CorrectRepetitionDORK(tmp,nav_ref,nav_range,params::Dict{Symbol,Any})
    
    if params[:kz_enc] == 0  # Linear
        nav = tmp[nav_range,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
        nav = mean(nav, dims=2)
    elseif params[:kz_enc] == 1 # Center-out
        nav = tmp[nav_range,1,:,1]
        nav = mean(nav, dims=2)
    end

    @info("Stop.. rDORK...")
    @infiltrate

    del_omg = mean(angle.(nav)-angle.(nav_ref))/params[:TE]

    tmp = tmp.*exp.(-1im.*del_omg.*params[:acq_times]')

    return tmp
end

"""
CorrectInterleafDORK(tmp,nav_range,Interleaf::Int,params::Dict{Symbol,Any})

Performs Interleaf DORK correction
"""

function CorrectInterleafDORK(tmp,nav_range,Interleaf::Int,params::Dict{Symbol,Any})

    if params_pulseq["gen"]["kz_enc"] == 0 # Linear
        # Here the reference interleaf is the first one
        nav_ref = tmp[nav_range,Int(numPar/2+1),:,1]
        nav_ref = mean(nav_ref, dims=2)

        nav = tmp[nav_range.*Interleaf.+numRead,Int(numPar/2+1),:,1]
        nav = mean(nav, dims=2)
    else params_pulseq["gen"]["kz_enc"] == 1 # Center-out
        # Here the reference interleaf is the first one
        nav_ref = tmp[nav_range,1,:,1]
        nav_ref = mean(nav_ref, dims=2)

        nav = tmp[nav_range.*Interleaf.+numRead,1,:,1]
        nav = mean(nav, dims=2)
    end

    del_omg = mean(angle.(nav)-angle.(nav_ref))/params[:TE]

    tmp[numRead*Interleaf+1:numRead*Interleaf+numRead,:,:,:] = tmp[numRead+1:numRead+numRead,:,:,:].*exp.(-1im.*del_omg.*params[:acq_times][1:numRead])

    return tmp
end


"""
Correctk0(RawData::RawAcquisitionData,k0::AbstractVector,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function Correctk0(tmp,k0_meas::AbstractMatrix{Float64},k0_sim::AbstractMatrix{Float64},params::Dict{Symbol,Any})
    
    # Un-do Siemens ECC
    tmp = tmp./exp.(1im.*k0_sim)

    # Apply Skope/girf K0
    tmp = tmp.*exp.(1im.*k0_meas.*-2π)

    return tmp
end


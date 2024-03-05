
"""
CorrectPartitionDORK(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function CorrectPartitionDORK(tmp,nav_ref,params::Dict{Symbol,Any})
    
    if params[:kz_enc] == 0  # Linear
        nav = tmp[6:30,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
        nav = mean(nav, dims=2)
    elseif params[:kz_enc] == 1 # Center-out
        nav = tmp[6:30,1,:,1]
        nav = mean(nav, dims=2)
    end

    del_phi = mean(angle.(nav)-angle.(nav_ref))/params[:TE]
    
    tmp = tmp.*exp.(-1im.*del_phi.*params[:acq_times]')

    return tmp
end


"""
Correctk0(RawData::RawAcquisitionData,k0::AbstractVector,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function Correctk0(tmp,k0::AbstractVector{Float64},params::Dict{Symbol,Any})
    
    k0 = reshape(k0,size(tmp)[1],size(tmp)[2])
    tmp = tmp.*exp.(1im.*-k0)

    return tmp
end
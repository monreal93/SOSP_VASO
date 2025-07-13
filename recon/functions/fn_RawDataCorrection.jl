
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

Performs repetition DORK correction
"""

function CorrectRepetitionDORK(tmp,nav,nav_ref,params::Dict{Symbol,Any})

    @info("Stop... Repetition DORK...")
    @infiltrate

    if nav === nothing
        if params[:kz_enc] == 0  # Linear
            nav = tmp[:,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
            nav = mean(nav, dims=2)
        elseif params[:kz_enc] == 1 # Center-out
            nav = tmp[:,1,:,1]
            nav = mean(nav, dims=2)
        end
    else
        # nav = nav[:,1,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1,1]  # Only center partiton and first FID
        # nav_ref = nav_ref[:,1,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:]  # Only center partition
        nav = nav[:,1,:,:,1,1]
    end

    del_omg = mean(DSP.unwrap(angle.(nav), dims=1)-DSP.unwrap(angle.(nav_ref),dims=1),dims=1)/params[:TE]

    tmp = tmp.*exp.(-1im.*del_omg.*params[:acq_times]')

    return tmp
end

"""
CorrectInterleafDORK(tmp,nav_range,Interleaf::Int,params::Dict{Symbol,Any})

Performs Interleaf DORK correction
"""

function CorrectInterleafDORK(tmp,nav,params::Dict{Symbol,Any})

    @info("Stop... Interleaf DORK...")
    @infiltrate

    # Get interl dim to the end
    nav = permutedims(nav,[1,3,4,2,5])
    
    del_phi_n_N = DSP.unwrap(mean(DSP.unwrap(angle.(nav[:,:,1,1,:]),dims=1).-DSP.unwrap(angle.(nav[:,:,1,1,1] ),dims=1), dims=1), dims=2) 
    # del_phi_n_I = DSP.unwrap(mean(DSP.unwrap(angle.(nav[:,:,1,2,:]),dims=1).-DSP.unwrap(angle.(nav[:,:,1,2,1] ),dims=1), dims=1), dims=2)

    # omg_n_I = (del_phi_n_I./ params_pulseq["gen"]["acqTR"])
    omg_n_N = (del_phi_n_N./ params_pulseq["gen"]["acqTR"])
    # del_omg_n = (omg_n_I.-omg_n_N)
    del_omg_n = omg_n_N

    # Get interl dim to second dim
    del_omg_n = permutedims(del_omg_n,[1,3,2])
    del_omg_n = repeat(del_omg_n,Int.(params_pulseq["gen"]["ro_samples"]),1,1)
    del_omg_n = reshape(del_omg_n,size(tmp,1),size(del_omg_n,3))

    tmp = tmp.*exp.(-1im.*del_omg_n.*params[:acq_times]')

    return tmp
end

"""
CorrectPartitionDORK(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function CorrectPartitionDORK(tmp,nav,params::Dict{Symbol,Any})

    @info("Stop... Partition DORK...")
    @infiltrate

    # del_phi_n = mean(DSP.unwrap(angle.(nav[:,1,:,1,2]),dims=1).-DSP.unwrap(angle.(nav[:,1,:,1,1]),dims=1), dims=1) ./-params[:fid_nav_ΔTE]
    # del_phi_r = del_phi_n[:,1]
    # del_omg = (del_phi_n.-del_phi_r)

    # del_phi_n_N = DSP.unwrap(mean(DSP.unwrap(angle.(nav[:,1,:,1,1]),dims=1).-DSP.unwrap(angle.(nav[:,1,1,1,1] ),dims=1), dims=1), dims=2) 
    # del_phi_n_I = DSP.unwrap(mean(DSP.unwrap(angle.(nav[:,1,:,1,2]),dims=1).-DSP.unwrap(angle.(nav[:,1,1,1,2] ),dims=1), dims=1), dims=2)

    # omg_n_I = (del_phi_n_I./ params_pulseq["gen"]["acqTR"]*params_pulseq["spi"]["interl"])
    # omg_n_N = (del_phi_n_N./ params_pulseq["gen"]["acqTR"]*params_pulseq["spi"]["interl"])
    # del_omg_n = (omg_n_I.-omg_n_N) # ./ params[:fid_nav_ΔTE] 
    # del_phi_n_o = ((2e-3*del_phi_n_N) -  (2.7e-3*del_phi_n_I)) ./ params[:fid_nav_ΔTE]
    # del_omg_n = omg_n_N

    #### First nav, then partitions
    # nav = mean(nav, dims=4)   # Mean over channels..
    del_omg_n = nav[:,1,:,:,1,2] .* conj.(nav[:,1,:,:,1,1])
    # del_omg_n = nav[:,1,:,:,1,2] ./ (nav[:,1,:,:,1,1])
    del_omg_n = DSP.unwrap(angle.(del_omg_n),dims=1) ./ params[:fid_nav_ΔTE]
    del_omg_n = mean(del_omg_n, dims=1)   # Mean over readout
    # del_omg_n = mean(del_omg_n, dims=3)   # Mean over channels..
    del_omg_n = permutedims(del_omg_n,(3,2,1))
    del_omg_n = del_omg_n[:,:,1] .- del_omg_n[:,1,1]  # Reference partition
    del_omg_n = reshape(del_omg_n, (size(del_omg_n)...,1))
    del_omg_n = permutedims(del_omg_n,(3,2,1))

    # Temp.. setting to zero a channel that is bad..
    # del_omg_n[:,:,13] = del_omg_n[:,:,13].*0 

    # #### First partitions, then nav
    # del_phi_n = DSP.unwrap(angle.(nav[:,1,:,:,1,:]), dims=1)
    # del_phi_R = DSP.unwrap(angle.(nav[:,1,1,:,1,:]), dims=1)  # Reference partition..
    # del_phi_R = reshape(del_phi_R,size(del_phi_R)[1]... ,1... ,size(del_phi_R)[2]...,size(del_phi_R)[3])
    # del_phi_n_N = del_phi_n[:,:,:,1] .- (del_phi_R[:,:,:,1])
    # del_phi_n_I = del_phi_n[:,:,:,2] .- (del_phi_R[:,:,:,2])
    # del_omg_n = (del_phi_n_I .- del_phi_n_N)./ params[:fid_nav_ΔTE]

    # del_omg = (del_phi_i-del_phi_n)./params[:fid_nav_ΔTE]
    # del_phi_o = (-params[:fid_nav_ΔTE]*del_phi_i)./-params[:fid_nav_ΔTE]
    # del_phi_o = del_phi_o .* (-1)
    # del_t = maximum(params[:acq_times]) - minimum(params[:acq_times])
    # del_omg = (del_phi_i-del_phi_n) ./ del_t

    # tmp = tmp.*exp.(-1im.*(del_phi_n_o.+del_omg_n).*params[:acq_times]')
    tmp = tmp.*exp.(1im.*del_omg_n.*params[:acq_times]')

    return tmp
end


"""
Correctk0(RawData::RawAcquisitionData,k0::AbstractVector,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function Correctk0(tmp,k0_meas::AbstractMatrix{Float64},k0_sim::Any,params::Dict{Symbol,Any})

    # Un-do Siemens ECC
    if params[:traj_type] == "sk"
        # tmp = tmp./exp.(1im.*k0_sim.*2*pi)
        tmp = tmp./exp.(1im.*k0_sim)
    end

    # Apply Skope/girf K0
    tmp = tmp.*exp.(1im.*k0_meas.*2π)
    # tmp = tmp.*exp.(-1im.*k0_meas.*1)

    return tmp
end


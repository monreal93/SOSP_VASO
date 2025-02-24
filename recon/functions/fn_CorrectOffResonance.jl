using FLoops
using NearestNeighbors

"""
CorrectOffResonanceTimeSegmented(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

"""
function CorrectOffResonanceTimeSegmented(ks_traj,kdata,OffResonanceMap,params,params_recon)

    # Get OffResonanceMap
    OffResonanceMap = imag.(OffResonanceMap)./2π
    # OffResonanceMap = OffResonanceMap .- maximum(OffResonanceMap)  # AMM: Should I have positive values?? 

    # Temp: trying stuff
    ks_traj.times = ks_traj.times .- ks_traj.times[1]

    # Time Segments
    Times = ks_traj.times[1:Int(params["gen"]["ro_samples"])]
    TimeSegments = Int(round(4*(maximum(OffResonanceMap)-minimum(OffResonanceMap))*params["gen"]["ro_time"]))
    TimeSegments = Int(ceil(TimeSegments/2)*2)
    # Temp: Increasing the number of time segments:
    # TimeSegments += 10

    SegmentLength = round(params["gen"]["ro_samples"]/TimeSegments/2)*2

    # Empty Ireco
    Ireco = Array{ComplexF32}(undef,Tuple([params["gen"]["n"]... TimeSegments]))

    @time @floop for i_TimeSegments in 1:TimeSegments
    # for i_TimeSegments in 1:Int(TimeSegments)
        if i_TimeSegments == TimeSegments
            SegmentRange = Int.(SegmentLength*(i_TimeSegments-1)+1:params["gen"]["ro_samples"])
        else
            SegmentRange = Int.(SegmentLength*(i_TimeSegments-1)+1:SegmentLength*(i_TimeSegments))
        end
        
        SegmentTimes = Times[SegmentRange]
        SegmentIndex = findall(x -> x in SegmentTimes, ks_traj.times)

        # Spliting trajectory and data into segments
        kdata_seg = deepcopy(kdata)
        kdata_seg[1] = kdata_seg[1][SegmentIndex,:]
        ks_traj_seg = deepcopy(ks_traj)
        ks_traj_seg.nodes = ks_traj_seg.nodes[:,SegmentIndex]
        ks_traj_seg.times = ks_traj_seg.times[SegmentIndex]
        ks_traj_seg.numSamplingPerProfile = Int(SegmentLength)

        # Create AcqData objects
        acqData = AcquisitionData(ks_traj_seg,kdata_seg; encodingSize=(Int.(SegmentLength),1,Int(params["gen"]["n_ov"][3])), fov=Tuple(params["gen"]["fov"]))

        # Reconstruction
        # @info(string("Reconstructing Repetition # ",i_rep))
        Ireco_seg = reconstruction(acqData, params_recon)

        # @info("Stop... Time-Segmented B0...")
        # @infiltrate

        local Ireco[:,:,:,i_TimeSegments] = Ireco_seg.data[:,:,:,1,1,1] .* exp.(im*2π.*OffResonanceMap.*SegmentTimes[1])

    end

    Ireco = dropdims(sum(Ireco,dims=4),dims=4)

    # @info("Stop... Time-Segmented B0...")
    # @infiltrate

    return Ireco

end

"""
CorrectOffResonanceFrequencySegmented(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

"""
function CorrectOffResonanceFrequencySegmented(ks_traj,kdata,OffResonanceMap,params,params_recon)

    # ks_traj.times = ks_traj.times .- ks_traj.times[1]
    Times = ks_traj.times

    # Get OffResonanceMap
    OffResonanceMap = imag.(OffResonanceMap)./2π
    # OffResonanceMap = OffResonanceMap .- maximum(OffResonanceMap)  # AMM: Should I have positive values?? 
    # OffResonanceMap = OffResonanceMap .* (2π)

    OffResonanceMapRange = zeros(size(OffResonanceMap))

    # Frequency Segments
    FrequencySegments = Int(round(4*(maximum(OffResonanceMap)-minimum(OffResonanceMap))*params["gen"]["ro_time"]))
    FrequencySegments = Int(ceil(FrequencySegments/2)*2)
    # FrequencySegments = 40
    SegmentRanges = LinRange(minimum(OffResonanceMap),maximum(OffResonanceMap),FrequencySegments)

    for i_FrequencySegments in 1:Int(FrequencySegments)-1
        SegmentIndex = findall(x -> SegmentRanges[i_FrequencySegments] <= x <= SegmentRanges[i_FrequencySegments+1], OffResonanceMap)
        OffResonanceMapRange[SegmentIndex] .= (SegmentRanges[i_FrequencySegments]+SegmentRanges[i_FrequencySegments+1])/2
    end

    SegmentValues = sort(unique(OffResonanceMapRange))

    # Empty Ireco
    Ireco = zeros(ComplexF32, size(OffResonanceMap))

    # @info("Stop... Frequency-Segmented B0...")
    # @infiltrate

    @time @floop for i_FrequencySegments in 1:Int(FrequencySegments)-1
    # for i_FrequencySegments in 1:Int(FrequencySegments)-1

        # Adding phase to segment
        kdata_seg = deepcopy(kdata)
        kdata_seg[1] = kdata_seg[1] .* exp.(im*2π.*SegmentValues[i_FrequencySegments].*Times)

        # Create AcqData objects
        acqData = AcquisitionData(ks_traj,kdata_seg; encodingSize=(Int(params["gen"]["ro_samples"]),1,Int(params["gen"]["n_ov"][3])), fov=Tuple(params["gen"]["fov"]))

        # Reconstruction
        Ireco_seg = reconstruction(acqData, params_recon)
        Ireco_seg = dropdims(Ireco_seg,dims=(4,5,6))

        SegmentIndex = OffResonanceMapRange .== SegmentValues[i_FrequencySegments]

        local Ireco[SegmentIndex] = Ireco_seg.data[SegmentIndex]

        # @info("Stop... Frequency-Segmented B0...")
        # @infiltrate

        # Trying to clear some memory
        ccall(:malloc_trim, Cvoid, (Cint,), 0) 
        GC.gc()

    end

    return Ireco

end
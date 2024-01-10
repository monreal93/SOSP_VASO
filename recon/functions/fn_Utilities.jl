"""
MergeReconstrucion(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Merge each repetition of the fMRI series into a timeseries
"""
function MergeReconstrucion(path::String,scan::String,Repetitions::UnitRange,ParallelRecon::Bool,OffResonanceCorrection::Bool)

    file_prefix = string(path,"/recon/3d/",scan)
    file_suffix = ""
    if ParallelRecon
        file_suffix = string(file_suffix,"_cs")
    end
    if OffResonanceCorrection
        file_suffix = string(file_suffix,"_b0")
    end

    # Read one, to get the volume dimensions
    tmp = niread(string(file_prefix,"_rep_",Repetitions[1],file_suffix,".nii"))

    @infiltrate
    TimeSeries = Array{Float32,4}(undef,(size(tmp.raw)[1:3]...,length(Repetitions)...))
    for i_rep=1:length(Repetitions)
        tmp = niread(string(file_prefix,"_rep_",Repetitions[i_rep],file_suffix,".nii"))
        TimeSeries[:,:,:,i_rep] .= tmp.raw[:,:,:,1]
    end

    # Scailing the data, for now arbitrary
    TimeSeries = TimeSeries.*1e9 

    return TimeSeries
end
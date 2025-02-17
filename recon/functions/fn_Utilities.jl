"""
MergeReconstrucion(path::String,scan::String,Repetitions::Any,ParallelRecon::Bool,OffResonanceCorrection::Bool,CocoCorrection::Bool, k0Correction::Bool,rDORKCorrection::Bool,hOrderCorrection::Bool; TrajectoryType::String="", Suffix="",B0CorrectionType)

Merge each repetition of the fMRI series into a timeseries
"""
function MergeReconstrucion(path::String,scan::String,Repetitions::Any,ParallelRecon::Bool,OffResonanceCorrection::Bool,CocoCorrection::Bool, k0Correction::Bool,rDORKCorrection::Bool,hOrderCorrection::Bool; TrajectoryType::String="nom", Suffix="",B0CorrectionType)

    file_prefix = string(path,"/recon/3d/",scan)
    file_suffix = string("_",TrajectoryType) 

    if ParallelRecon
        file_suffix = string(file_suffix,"_cs")
    end
    if OffResonanceCorrection
        if hOrderCorrection
            file_suffix = string(file_suffix,"_hb0")
        end
        if B0CorrectionType==0
            file_suffix = string(file_suffix,"_b0")
        elseif B0CorrectionType==1
            file_suffix = string(file_suffix,"_tsb0")
        elseif B0CorrectionType==2
            file_suffix = string(file_suffix,"_fsb0")
        end
    end
    if CocoCorrection
        file_suffix = string(file_suffix,"_co") 
    end
    if k0Correction
        file_suffix = string(file_suffix,"_k0")
    end
    if rDORKCorrection
        file_suffix = string(file_suffix,"_rDORK")
    end
    file_suffix = string(file_suffix,Suffix)

    # Read one, to get the volume dimensions
    tmp = niread(string(file_prefix,"_rep_",Repetitions[1],file_suffix,".nii"))

    TimeSeries = Array{Float32,4}(undef,(size(tmp.raw)[1:3]...,length(Repetitions)...))
    for i_rep=1:length(Repetitions)
        tmp = niread(string(file_prefix,"_rep_",Repetitions[i_rep],file_suffix,".nii"))
        TimeSeries[:,:,:,i_rep] .= tmp.raw[:,:,:,1]
    end

    # Scailing the data, for now arbitrary
    TimeSeries = TimeSeries.*1e11 

    return TimeSeries
end

"""
MergeReconstrucion(path::String,scan::String,suffix::String,Repetitions::Any)

Merge each repetition of the fMRI series into a timeseries
"""
function MergeReconstrucion(path::String,scan::String,suffix::String,Repetitions::Any)

    file_prefix = string(path,"/recon/3d/",scan)

    # Read one, to get the volume dimensions
    tmp = niread(string(file_prefix,"_rep_",Repetitions[1],suffix,".nii"))

    TimeSeries = Array{Float32,4}(undef,(size(tmp.raw)[1:3]...,length(Repetitions)...))
    for i_rep=1:length(Repetitions)
        tmp = niread(string(file_prefix,"_rep_",Repetitions[i_rep],suffix,".nii"))
        TimeSeries[:,:,:,i_rep] .= tmp.raw[:,:,:,1]
    end

    # Scailing the data, for now arbitrary
    TimeSeries = TimeSeries.*1e11 

    return TimeSeries
end
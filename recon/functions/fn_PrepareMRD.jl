

"""
SetRepetitionLabel(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Adds repetition label to raw data object
"""
function SetRepetitionLabel(rawData,params)
    numPar = params[:numPar]
    numSet = params[:numSet]
    
    # Here I want to loop over the profiles and add the "repeition" labels, need to loop over parittions... and adc split 
    i_rep = UInt16(0)
    for i_profile=1:length(rawData.profiles)
        rawData.profiles[i_profile].head.idx.repetition = i_rep
        rawData.profiles[i_profile].head.idx.slice = rawData.profiles[i_profile].head.idx.kspace_encode_step_2
        # Repetitions
        if mod(i_profile,numPar*numSet) == 0
            i_rep += UInt16(1)
        end
    end

    return rawData
end

"""
    FormatRawData(RawData::RawAcquisitionData;single_rep::bool=false,rep_format::Int))

Rearanges the raw data, into [readouts*set,partitions,channels,repetitions]
"""
function FormatRawData(rawData,params;single_rep::Bool=false,rep_format::Int,fid_nav::Int=0,MultiEchoGRE::Bool=false)
    numRead = params[:numRead]
    numInterl = params[:numInterl]
    numPar = params[:numPar]
    numSet = params[:numSet]
    numCha = params[:numCha] 
    numRep = params[:numRep]
    numEch = params[:numEch]

    # ToDo: Add functionality for multi-interleaves
    if single_rep == true
        if params_pulseq["gen"]["seq"] == 3 && MultiEchoGRE == false
            # if mod(rep_format,numEch+1) == 0
            #     @info("Stop... Format RawData XXXXXX....")
            #     @infiltrate
            #     i_prof = Int(((rep_format-1)*numEch*numPar*numInterl*numSet)+1)
            # else
            #     i_prof = Int(((rep_format-1)*(numEch-1)*numInterl*numSet)+1)
            # end
            i_prof = Int(((mod(rep_format-1,numEch))*(numEch-1)*numInterl*numSet)+1) + Int(numEch*numPar*numInterl*numSet*numEch*(ceil((rep_format-numEch)/numEch)))
        else
            i_prof = Int(((rep_format-1)*numPar*numInterl*numSet)+1)
        end
        numRep = 1
        tmp = Array{ComplexF32}(undef,Int(numRead/numSet),numSet,numInterl,numPar,numCha,1)
        kdata = Array{ComplexF32}(undef,numRead*numInterl,numPar,numCha,1)
    else
        i_prof = 1
        tmp = Array{ComplexF32}(undef,Int(numRead/numSet),numSet,numInterl,numPar,numCha,numRep)
        kdata = Array{ComplexF32}(undef,numRead*numInterl,numPar,numCha,numRep)
    end

    if fid_nav == 1
        if params_pulseq["gen"]["seq"] == 3 && MultiEchoGRE == false
            i_prof = i_prof+(2*(mod(rep_format-1,numEch))*numSet) 
        else
            # i_prof = ((i_prof-1)*2)+1  # Original
            i_prof +=  ((rep_format-1)*numPar*numInterl*2)
        end
        ndata = Array{ComplexF32}(undef,size(rawData.profiles[1].data,1),numInterl,numPar,numCha,numRep,2)
    end

    for i_rep=1:numRep
        for i_par=1:numPar
            for i_interl=1:numInterl
                if fid_nav == 1
                    for i_nav=1:2
                        ndata[:,i_interl,i_par,:,i_rep,i_nav] = rawData.profiles[i_prof].data # Original
                        i_prof += 1
                    end
                end
                for i_set=1:numSet
                    tmp[:,i_set,i_interl,i_par,:,i_rep] = rawData.profiles[i_prof].data # Original
                    # tmp[:,i_set,i_interl,i_par,:,i_rep] = rawData.profiles[i_prof].data[1:numRead,:]
                    i_prof += 1
                end
            end
            if params_pulseq["gen"]["seq"] == 3 && MultiEchoGRE == false
                i_prof += (numInterl*(numEch))
                if fid_nav == 1
                    i_prof += (2*(numEch))
                end
            end
        end
    end

    kdata .= reshape(tmp,(:,numPar,numCha,numRep))

    if @isdefined ndata
        return kdata, ndata
    else
        return kdata, nothing
    end
end

"""
MergeTrajectory(rawData::RawAcquisitionData,ks_traj::....;skope::Bool=false)

Merges trajectory with 
"""
function MergeTrajectory(rawData,params)
    numSet = params[:set]
    numPar = params[:sl]

    # Here I want to loop over the profiles and add the "repeition" labels, need to loop over parittions... and adc split 
    i_rep = UInt16(0)
    i_par = UInt16(0)
    for i_profile=1:length(rawData.profiles)
        rawData.profiles[i_profile].head.idx.repetition = i_rep
        rawData.profiles[i_profile].head.idx.slice = rawData.profiles[i_profile].head.idx.kspace_encode_step_2
        # Repetitions
        if mod(i_profile,numPar*numSet) == 0
            global i_rep += UInt16(1)
        end
    end

    return rawData
end

"""
ReconCartesianData(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Performs cartesian reconstruction, FFT, either 2D or 3D
"""
function ReconCartesianData(acqData::AcquisitionData,dims::Int; interleaved=false)
    ## AMM: ToDo:
    # - Add functionality to select if FFT 2D or 3D...
    #

    numEchos = size(acqData.kdata,1)
    # numSlices = size(acqData.kdata,2)             # Slices or Phase Enc 2 in 3D
    numRead = acqData.encodingSize[1]
    numPart = acqData.encodingSize[2]
    numSlices = acqData.encodingSize[3]
    numChan = size(acqData.kdata[1],2)

    recon = Array{ComplexF32,5}(undef,numRead,numPart,numSlices,numChan,numEchos)

    if dims==2
        # Reshape the data and do FFT
        for i_ech=1:numEchos
            for i_sl=1:numSlices

                tmp = acqData.kdata[i_ech,i_sl,1]
                tmp = reshape(tmp,numRead,numPart,numChan)
                
                # # FFT, only for 2D now, FFTW use:
                # tmp = fftshift(ifft(tmp),1)
                # tmp = fftshift(ifft(tmp),2)

                # FFT, only for 2D now, MKL use:
                tmp = fftshift(ifft(tmp))
                # tmp = sqrt.(sum((abs.(tmp).^2),dims=3))

                recon[:,:,i_sl,:,i_ech] = tmp
            end
        end

    elseif dims==3
        for i_ech=1:numEchos
            tmp = acqData.kdata[i_ech,1,1]
            tmp = reshape(tmp,numRead,numSlices,numPart,numChan)
            
            # # FFT, only for 2D now, FFTW use:
            # tmp = fftshift(ifft(tmp),1)
            # tmp = fftshift(ifft(tmp),2)

            # FFT, only for 2D now, MKL use:
            tmp = fftshift(ifft(tmp))
            # tmp = sqrt.(sum((abs.(tmp).^2),dims=3))

            tmp = permutedims(tmp,[1,3,2,4])

            recon[:,:,:,:,i_ech] = tmp
        end
    end

    # Removing oversampling
    recon = recon[Int((end/2)-(size(recon,1)/4)):Int((end/2)+(size(recon,1)/4)-1),:,:,:,:]

    # Rearranging if needed.. only when interleaved acquisition...
    if interleaved
        tmp = Array{ComplexF32,5}(undef,size(recon))
        tmp[:,:,2:2:end,:,:] = recon[:,:,1:Int(size(recon,3)/2),:,:]
        tmp[:,:,1:2:end,:,:] = recon[:,:,Int(size(recon,3)/2)+1:end,:,:]
        recon = tmp
    end

    recon = reverse(recon, dims=1)

    return recon
end

"""
ReconCartesianData(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Performs cartesian reconstruction, FFT, either 2D or 3D, taking an Array as input
"""
function ReconCartesianData(rawData::Array{ComplexF64}; dims::Int=2)

    recon = Array{ComplexF64,}(undef,size(rawData))

    recon = fftshift(fft(rawData,1),1)
    recon = fftshift(fft(recon,2),2)

    if dims==3
        recon = fftshift(fft(recon,3),3)
    end

    # Convert to Float32
    recon = convert(Array{ComplexF32}, recon)

    return recon
end
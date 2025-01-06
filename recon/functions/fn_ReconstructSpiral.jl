
"""
ReconSpiralMEGRE(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Performs ME-GRE spiral reconstruction, simple regridding with NFFT
"""
function ReconSpiralMEGRE(acqData::AcquisitionData{Float32,3},params::Dict{String,Any})

    # params_pulseq["gen"]["n_ov"] = Int.(params["gen"]["n_ov"])

    params_recon = Dict{Symbol, Any}()
    params_recon = merge(defaultRecoParams(), params_recon)
    params_recon[:reconSize] = Int.((params["gen"]["n_ov"][1],params["gen"]["n_ov"][2],params["gen"]["n_ov"][3]))
    params_recon[:regularization] = "L1"
    params_recon[:λ] = 1.e-2  # 1.e-2
    params_recon[:iterations] = 40 # (10)
    params_recon[:solver] = "admm" # cgnr (L2), fista (L1), admm(L1)
    params_recon[:method] = "nfft" # nfft, leastsquare

    recon = Array{ComplexF64}(undef,Int.(params_recon[:reconSize]))
    # Convert to Float32
    recon = convert(Array{ComplexF32,3}, recon)

    recon = reconstruction(acqData, params_recon)

    recon = recon[:,:,:,1,:,1]

    # Coil combine
    recon_sos = sqrt.(((sum(abs.(recon).^2; dims=4))))

    return recon
end

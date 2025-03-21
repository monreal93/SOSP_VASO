using MRIFieldmaps: b0map, b0model, b0init, b0scale
using ROMEO

"""
CalculateSensitivitityMap(RawData::RawAcquisitionData,params::Dict{Symbol,Any})
MtxSize = Tuple (nx,ny,nz) = Final size of Map
It expects raw data to be 3D
Calculates Sensitivity map
"""
function CalculateSensitivityMap(recon,MtxSize::Tuple;calib_size::Int=12)

    SensitivityMap = Array{ComplexF32}(undef,MtxSize)
    calibration = Array{ComplexF32}(undef,size(recon)[1:4])

    numChan = size(recon)[4]

    for i_ch = 1:numChan
        calibration[:,:,:,i_ch] = fftshift(fftshift(fftshift(fft(recon[:,:,:,i_ch,1]),1),2),3)
    end

    # # Back to k-space
    # calibration = recon[:,:,:,:,1]
    # calibration = fft(calibration,1)
    # calibration = fft(calibration,2)
    # calibration = fftshift(fft(calibration,3),3)

    if MtxSize[3] > calib_size*2
        calibration = calibration[Int(round(end/2+1-calib_size)):Int(round(end/2+calib_size)),Int(round(end/2+1-calib_size)):Int(round(end/2+calib_size)),Int(floor(end/2+1)-calib_size):Int(floor(end/2)+calib_size),:]
    else
        calibration = calibration[Int(round(end/2+1-calib_size)):Int(round(end/2+calib_size)),Int(round(end/2+1-calib_size)):Int(round(end/2+calib_size)),:,:]
    end
        ## ESPIRIT wants even numbers in slice direction, lets add some slices
    # if mod(Int(MtxSize[3]./2),2) == 1
    #     MtxSize_tmp = (MtxSize[1],MtxSize[2],Int(MtxSize[3]+2))
    # else
    #     MtxSize_tmp = copy(MtxSize)
    # end

    SensitivityMap = espirit(calibration,MtxSize)
    SensitivityMap = fftshift(fftshift(SensitivityMap,1),2)

    # ToDo: Sometimes I need to do fftshift in 3 dim also
    SensitivityMap = fftshift(SensitivityMap,3)

    SensitivityMap = dropdims(SensitivityMap; dims = 5)

    # # ESPIRIT wants even numbers in slice direction, lets remove the added slices
    # SensitivityMap = SensitivityMap[:,:,2:end-1,:]

    return SensitivityMap
    
end

"""
CalculateOffresonanceMap(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Calculates B0 map
"""
function CalculateOffresonanceMap(recon_b0,SensitivityMap,EchoTimes::Vector{Float64})

    mask = SensitivityMap[:,:,:,1]
    mask[abs.(mask) .> 0] .= 1
    mask = imresize(mask,(size(recon_b0)[1],size(recon_b0)[2],size(recon_b0)[3]))
    mask = isone.(mask)

    # Lowering resolution of sensitivities
    smap = imresize(SensitivityMap,size(recon_b0)[1:4])

    # Taking only the first 3 echos
    ydata = recon_b0[:,:,:,:,1:3] .*1e5
    echotime = EchoTimes[1:3].*1e-3 *1s   # Original

    # ToDo: Do I need to flip diemension?
    # ydata = reverse(ydata, dims=1)

    yik_sos = sum(conj(smap) .* ydata; dims=4) # coil combine
    yik_sos = yik_sos[:,:,:,1,:] # (dims..., ne)

    (yik_sos_scaled, scale) = b0scale(yik_sos, echotime) # todo
    yik_scale = ydata / scale # original

    finit = b0init(ydata, echotime; smap)

    # ## Trying ROMEO phase unwrapping... Naive approach (1)
    # ph = angle.(yik_sos)
    # mag = abs.(yik_sos)
    # te = ustrip.(echotime) 
    # ph_unwrapped = ROMEO.unwrap(ph; mag=mag, TEs=te)
    # b0_romeo1 = (ph_unwrapped[:,:,:,1].-ph_unwrapped[:,:,:,2])./(((te[2]-te[1])))
    # b0_romeo2 = (ph_unwrapped[:,:,:,2].-ph_unwrapped[:,:,:,3])./(((te[3]-te[2])))
    # b0_romeo = (b0_romeo1.+b0_romeo2)./2
    # finit = b0_romeo .* s^-1

    ## Trying ROMEO phase unwrapping..(2)
    mag = abs.(yik_sos)
    te = ustrip.(echotime) 
    del_phi = angle.(recon_b0[:,:,:,:,1] .* conj.(recon_b0[:,:,:,:,2]))
    del_phi = ROMEO.unwrap(del_phi; mag=mag[:,:,:,1:2], TEs=te[1:2])
    del_phi = sum(del_phi, dims=4)
    finit = del_phi./(((te[2]-te[1]).*1e3))
    finit = finit[:,:,:,1] .* s^-1

    # # Original "Low smoothing":
    # fmap_run = (niter, precon, track; kwargs...) ->
    #     b0map(yik_scale, echotime; finit, smap, mask,
    #     order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)

    @info ("Stop... B0 map...")
    @infiltrate

    # Different smoothing values:
    # Default order =1
    fmap_run = (niter, precon, track; kwargs...) ->
    b0map(yik_scale, echotime; finit, smap, mask,
    order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)

    # # A lot of smoothing:
    # fmap_run = (niter, precon, track; kwargs...) ->
    #     b0map(yik_scale, echotime; finit, smap, mask,
    #     order=1, l2b=4, gamma_type=:PR, niter, precon, track, kwargs...)


    function runner(niter, precon; kwargs...)
        (fmap, times, out) = fmap_run(niter, precon, true; kwargs...)
        return (fmap, out.fhats, out.costs, times)
    end;
    # if !@isdefined(fmap_cg_d)
        niter_cg_d = 100  # (400)
        (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag)
        # (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :ichol)
    # end

    b0 = fmap_cg_d

    b0 = b0.*2π.*-1    # Original 
    # b0 = b0.*-1
    # b0[b0.≠0] = b0[b0.≠0] .-(20*2π)*1s^-1  # Temp: only for sv_01
    
    b0 = ustrip.(b0)

    b0 = imresize(b0,size(SensitivityMap)[1:3],method=BSpline(Cubic()))
    b0 = 1im.*b0
    b0 = convert(Array{ComplexF32,3},b0)

    @info ("Stop... B0 map...")
    @infiltrate

    return b0
end



"""
CalculateConcomitantFieldMap(RotMtx::AbstractArray{Float32},params::Dict{String,Any})

Calculates concomitant field map
"""
function CalculateConcomitantFieldMap(RotMatrix::Matrix{Float32},CenterPosition::Tuple,params_pulseq::Dict{String,Any})
    
    b0 = params_pulseq["gen"]["field_strength"]
    gm = params_pulseq["spi"]["max_grad"] .* 1e-3
    x = range(params_pulseq["gen"]["fov"][1]/2*-1,params_pulseq["gen"]["fov"][1]/2,Int(params_pulseq["gen"]["n_ov"][1]))
    y = range(params_pulseq["gen"]["fov"][2]/2*-1,params_pulseq["gen"]["fov"][2]/2,Int(params_pulseq["gen"]["n_ov"][2]))
    z = range(params_pulseq["gen"]["fov"][3]/2*-1,params_pulseq["gen"]["fov"][3]/2,Int(params_pulseq["gen"]["n_ov"][3]*params_pulseq["gen"]["kz"]))

    x = x .+ CenterPosition[1] 
    y = y .+ CenterPosition[2] 
    z = z .+ CenterPosition[3] 

    a1=RotMatrix[1,1]
    a2=RotMatrix[1,2]
    a3=RotMatrix[1,3]
    a4=RotMatrix[2,1]
    a5=RotMatrix[2,2]
    a6=RotMatrix[2,3]
    a7=RotMatrix[3,1]
    a8=RotMatrix[3,2]
    a9=RotMatrix[3,3]

    F1 = ((1/4)*((a1^2)+(a4^2))*((a7^2)+(a8^2))) + ((a7^2)*((a2^2)+(a5^2))) - (a7*a8*((a1*a2)+(a4*a5)))
    F2 = ((1/4)*((a2^2)+(a5^2))*((a7^2)+(a8^2))) + ((a8^2)*((a1^2)+(a4^2))) - (a7*a8*((a1*a2)+(a4*a5)))
    F3 = ((1/4)*((a3^2)+(a6^2))*((a7^2)+(a8^2))) + ((a9^2)*((a1^2)+(a2^2)+(a4^2)+(a5^2))) - (a7*a9*((a1*a3)+(a4*a6))) - (a8*a9*((a2*a3)+(a5*a6)))
    F4 = ((1/2)*((a2*a3)+(a5*a6))*((a7^2)-(a8^2))) + ((a8*a9)*((2*(a1^2))+(a2^2)+((2*(a4^2))+(a5^2)))) - (a7*a8*((a1*a3)+(a4*a6))) - (a7*a9*((a1*a2)+(a4*a5)))
    F5 = ((1/2)*((a1*a3)+(a4*a6))*((a8^2)-(a7^2))) + ((a7*a9)*((a1^2)+(2*(a2^2))+((a4^2)+(2*(a5^2))))) - (a7*a8*((a2*a3)+(a5*a6))) - (a8*a9*((a1*a2)+(a4*a5)))
    F6 = ((-1/2)*((a1*a2)+(a4*a5))*((a7^2)+(a8^2))) + ((a7*a8)*((a1^2)+(a2^2)+(a4^2)+(a5^2)))
    
    fc = Array{Float64}(undef,length(x),length(y),length(z))
    if sum(RotMatrix) == 3
        # Concomitant field (TRA,SAG,COR)
        for i_x=eachindex(x)
            for i_y=eachindex(y)
                for i_z=eachindex(z)    
                    fc[i_x,i_y,i_z] =  42.58e6 * ((gm^2)/(4*b0)) * (((x[i_x].^2)/4)+((y[i_y].^2)/4)+((z[i_z].^2)))
                end
            end
        end
    else
        # Rotated concomitant field (tilted)
        for i_x=eachindex(x)
            for i_y=eachindex(y)
                for i_z=eachindex(z)    
                    fc[i_x,i_y,i_z] =  42.58e6 * ((gm^2)/(4*b0)) * ((F1*(x[i_x].^2)) + (F2*(y[i_y].^2)) + (F3*(z[i_z].^2)) + (F4.*y[i_y].*z[i_z]) + (F5.*x[i_x].*z[i_z]) + (F6.*x[i_x].*y[i_y]))
                end
            end
        end
    end

    # Convert to rad/s
    CocoFieldMap = fc .* 2π

    # @info("Stop... CocoFieldMap...")
    # @infiltrate

    return CocoFieldMap
end
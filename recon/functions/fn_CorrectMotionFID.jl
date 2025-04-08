using CoordinateTransformations, Rotations, Distributions

function getMotionCalibrationMatrix(recon_b0,SensitivityMap,nav_ref;calib_steps=500)

    n_channels = size(SensitivityMap,4)
    n_times = calib_steps

    # Coil combine
    recon = recon_b0[:,:,:,:,1]       # Taking one of the echos...
    cs = imresize(SensitivityMap,Tuple(size(recon_b0)[1:4]))
    recon_b0_cc = sum(conj.(cs) .* recon; dims=4)
    recon_b0_cc = dropdims(recon_b0_cc, dims=4)

    # Random translations, from -4 to 4 mm
    d = Uniform(-10,10)
    tra = rand(d,3,calib_steps) .*1e-3
    # Random roation angles, from -4 to 4 deg
    d = Uniform(-10,10)
    α = rand(d,3,calib_steps)
    α = α.*pi/180

    # motion_params = [traX traY traZ pitch roll yaw; calib_steps]
    motion_params = [tra ; α]
    # n_terms = calib_steps.*6

    # nav = Vector{ComplexF32}(undef,n_terms)
    fid = Matrix{Float32}(undef,n_channels,n_times)
    calib = Matrix{ComplexF32}(undef,n_channels,6)


    @info("Stop... Motion corr...")
    @infiltrate

    # Find scaling factor...
    fid_ref = (recon_b0_cc.*cs)
    fid_ref = abs.(dropdims(sum(fid_ref,dims=(1,2,3)),dims=(1,2,3)))
    nav_scaling = abs.(nav_ref) ./ abs.(fid_ref)

    @floop for i=1:calib_steps

        t = Translation(motion_params[1,i],motion_params[2,i],motion_params[3,i])
        # r_pitch = recenter(RotMatrix3([1,0,0,0, cos(motion_params[4,i]), sin(motion_params[4,i]) , 0, -sin(motion_params[4,i]), cos(motion_params[4,i])]), center(recon))
        # r_roll = recenter(RotMatrix3([cos(motion_params[5,i]), 0, -sin(motion_params[5,i]), 0, 1, 0, sin(motion_params[5,i]), 0, cos(motion_params[5,i])]), center(recon))
        # r_yaw = recenter(RotMatrix3([cos(motion_params[6,i]), sin(motion_params[6,i]), 0, -sin(motion_params[6,i]), cos(motion_params[6,i]), 0, 0, 0, 1]), center(recon))

        r_pitch = recenter(RotMatrix3([1,0,0,0, cos(motion_params[4,i]), sin(motion_params[4,i]) , 0, -sin(motion_params[4,i]), cos(motion_params[4,i])]),  Tuple(size(recon_b0_cc)./2))
        r_roll = recenter(RotMatrix3([cos(motion_params[5,i]), 0, -sin(motion_params[5,i]), 0, 1, 0, sin(motion_params[5,i]), 0, cos(motion_params[5,i])]),  Tuple(size(recon_b0_cc)./2))
        r_yaw = recenter(RotMatrix3([cos(motion_params[6,i]), sin(motion_params[6,i]), 0, -sin(motion_params[6,i]), cos(motion_params[6,i]), 0, 0, 0, 1]),  Tuple(size(recon_b0_cc)./2))

        r = r_pitch ∘ r_roll ∘ r_yaw
        tr = t ∘ r

        img_trans = warp(recon_b0_cc, tr)
        img_trans = img_trans.parent
        replace!(img_trans, NaN=>0)

        # Cropping img_trans
        crop = Int.(floor.((size(img_trans) .- size(recon_b0_cc))./2))
        # crop_tmp = mod.(size(img_trans),2)
        crop_tmp = (isinteger.((size(img_trans) .- size(recon_b0_cc))./2))
        crop_tmp = ([Int(!crop_tmp[1]),Int(!crop_tmp[2]),Int(!crop_tmp[3])])
        img_trans = img_trans[crop[1]+crop_tmp[1]+1:end-crop[1],crop[2]+crop_tmp[2]+1:end-crop[2],crop[3]+crop_tmp[3]+1:end-crop[3]]
        # img_trans = img_trans[crop[1]+1:end-crop[1],crop[2]+1:end-crop[2],crop[3]+1:end-crop[3]]

        img_trans = (img_trans.*cs)

        fid[:,i] = abs.(dropdims(sum(img_trans,dims=(1,2,3)),dims=(1,2,3)))
        fid[:,i] = fid[:,i] .* nav_scaling

        img_trans = []

    end

    @info("Stop... Motion corr...")
    @infiltrate

    calib = fid / motion_params

    return calib
end

function applyMotionCalibrationMatrix(rawData,ksTraj::Matrix{Float32},nav,calib)

    @info("Stop... Apply Motion...")
    @infiltrate

    # AMM: Temp...
    nav = abs.(mean(nav,dims=3))
    nav = abs.(mean(nav,dims=1))
    nav = nav[:]
    # nav = nav[:] .* 10
    # nav = nav[1,1,21,:]

    motion_params = nav' / calib'
    # motion_params = abs.(motion_params)

    # Temp: manually selecting values
    # motion_params = [-1.377 -0.113 -0.004 -0.065 0.547 0.0005] # op1
    # motion_params = [-0.113 -0.004 -1.377 0.547 0.0005 -0.065] # op2
    # motion_params = [0 0 0 0 0 0] # op3
    # motion_params = [-0.113 -1.377 -0.004 0.547 -0.065 0.0005] # op4
    # motion_params = [-0.113 -1.377 -0.004 -0.065 0.547 0.0005] # op5 BEst one...
    # motion_params = [-0.1063 0.5812 -0.4436 0.6845 0.0847 0.0334]  # rep=240, sb_001_... 01092025_sb_9T

    motion_params[1:3] =  motion_params[1:3]./(params_pulseq["gen"]["res"].*1e3)'
    motion_params[4:6] =  motion_params[4:6]*(pi/180)

    # # Try 1, I got this sample values from the motion correction step of post-processing
    # motion_params[1:3] = [-0.7284,-0.0572,0.1293]./(params_pulseq["gen"]["res"].*1e3)'
    # motion_params[4:6] = [-0.1767*(pi/180), 0.1369*(pi/180), -0.246*(pi/180)]

    RotationMtx = RotXYZ(motion_params[4],motion_params[5],motion_params[6])

    ksTraj = RotationMtx * ksTraj

    # ksTraj = ksTraj[:,1:(params[:numRead]*params[:numInterl])]'
    ksTraj = permutedims(ksTraj,[2,1])
    ksTraj = reshape(ksTraj,Tuple([size(rawData)[1:2] ... 3]))

    rawData = rawData.*exp.(1im*2*π.*ksTraj[:,:,1].*(motion_params[1]))  # kx shift
    rawData = rawData.*exp.(1im*2*π.*ksTraj[:,:,2].*(motion_params[2]))  # ky shift
    rawData = rawData.*exp.(1im*2*π.*ksTraj[:,:,3].*(motion_params[3]))  # kz shift

    ksTraj = reshape(ksTraj,Tuple([size(rawData,1)*size(rawData,2) ... 3]))'

    # Normalizing trajectory
    ksTraj[1,:] =  ksTraj[1,:]./(maximum([abs(minimum(ksTraj[1,:])),abs(maximum(ksTraj[1,:]))])*2)
    ksTraj[2,:] =  ksTraj[2,:]./(maximum([abs(minimum(ksTraj[2,:])),abs(maximum(ksTraj[2,:]))])*2)
    ksTraj[3,:] =  ksTraj[3,:]./(maximum([abs(minimum(ksTraj[3,:])),abs(maximum(ksTraj[3,:]))])*2)

    @info("Stop... Apply Motion...")
    @infiltrate

    return rawData, ksTraj
end





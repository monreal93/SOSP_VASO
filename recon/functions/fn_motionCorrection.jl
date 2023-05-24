function motionCorrection()
# function motionCorrection(params::Dict{Symbol,Any},fm,cs)
    img = niread("/usr/share/sosp_vaso/data/04252023_sv_abc/tmp/gre.nii")
    fm = img.raw

    cs = matread("/usr/share/sosp_vaso/data/04252023_sv_abc/acq/cs_sv_01_238_236_24.mat"); cs = cs["coil_sens"];

    # mtx_s = [Int64(params[:mtx_s][1]),Int64(params[:mtx_s][2]),Int64(params[:mtx_s][3])]

    # fm = imresize(fm,(mtx_s[1],mtx_s[2],mtx_s[3],params[:nCoils]))
    fm = fm[:,:,:,1,:]

    img1 = sum(conj(cs) .* fm; dims=4) # coil combine
    img1 = dropdims(img1.*1e-2; dims=4)

    α = LinRange(-0.5,0.5,2)
    α = α*pi/180
    tra = LinRange(-0.1e-3,0.1e-3,2)
    n_terms = (length(α)^3)*(length(tra)^3)

    nav = Vector{ComplexF32}(undef,n_terms)
    motion_params = Matrix{Float32}(undef,n_terms,6)
    calib = Matrix{ComplexF32}(undef,32,n_terms)

    # idx = 1

    @infiltrate

    for i_x=eachindex(tra), i_y=eachindex(tra), i_z=eachindex(tra), i_pitch=eachindex(α), i_roll=eachindex(α), i_yaw=eachindex(α)
    # @floop for i_x=eachindex(tra), i_y=eachindex(tra), i_z=eachindex(tra), i_pitch=eachindex(α), i_roll=eachindex(α), i_yaw=eachindex(α)

            idx = (i_x^2) + (i_y^2) + (i_z^2) + (i_pitch^2) + (i_roll^2) + (i_yaw^2) -5 
            idx1 = ((i_x-1)^2) + ((i_y-1)^2) + ((i_z-1)^2) + ((i_pitch-1)^2) + ((i_roll-1)^2) + (i_yaw)

            @infiltrate
            # t = Translation(0,0,0)
            # r_x = recenter(RotMatrix3([1,0,0,0, cos(α), sin(α) , 0, -sin(α), cos(α)]), center(img1))
            # r_y = recenter(RotMatrix3([cos(α), 0, -sin(α), 0, 1, 0, sin(α), 0, cos(α)]), center(img1))
            # r_z = recenter(RotMatrix3([cos(α), sin(α), 0, -sin(α), cos(α), 0, 0, 0, 1]), center(img1))

            t = Translation(tra[i_x],tra[i_y],tra[i_z])
            r_pitch = recenter(RotMatrix3([1,0,0,0, cos(α[i_pitch]), sin(α[i_pitch]) , 0, -sin(α[i_pitch]), cos(α[i_pitch])]), center(img1))
            r_roll = recenter(RotMatrix3([cos(α[i_roll]), 0, -sin(α[i_roll]), 0, 1, 0, sin(α[i_roll]), 0, cos(α[i_roll])]), center(img1))
            r_yaw = recenter(RotMatrix3([cos(α[i_yaw]), sin(α[i_yaw]), 0, -sin(α[i_yaw]), cos(α[i_yaw]), 0, 0, 0, 1]), center(img1))

            r = r_pitch ∘ r_roll ∘ r_yaw
            tr = t ∘ r

            img_trans = warp(img1, tr, indices_spatial(img1))
            replace!(img_trans, NaN=>0)

            img_trans = img_trans.*cs
            img_tmp = reshape(img_trans,(:,32))

            # (string("Done with n=",idx))

            # nav[idx] = sum(img_trans[:])
            # @infiltrate
            motion_params[idx,:] = [tra[i_x],tra[i_y],tra[i_z],α[i_pitch],α[i_roll],α[i_yaw]] 
            calib[:,idx] = sum(img_tmp; dims=1)'
            
            # @infiltrate

            # idx += 1

            println("n = ", idx , " in Thread ID: ", Threads.threadid())
        end
@infiltrate
calib = abs.(calib)

return calib
end





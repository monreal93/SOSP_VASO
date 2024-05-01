function fn_sv_recon(params_sv::Dict{Symbol,Any})
    # if params_sv[:fmri] == 1
    #     file_name = string(params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl]);
    # else
    #     file_name = string(params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl]);
    # end

    file_name = string(params_sv[:fieldmap],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl]);
    
    if params_sv[:do_pi_recon]
        if isfile(string(params_sv[:path],"acq/cs_",file_name,".mat"))
            @info ("Sensitivity map exists... Re-calculate? y/n")
            input = readline()
        else
            input = ""
        end
        if input == "n"
            @info ("Loading Sensitivity maps ...")
            cs = matread(string(params_sv[:path],"acq/cs_",file_name,".mat")); cs = cs["coil_sens"]

            # Loading fieldmap and saving in image domain
            fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:fieldmap],"_fieldmap.mat"))

            # fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:scan],"_fieldmap.mat"))
            fieldmap_ft = fieldmap_ft["b0"]
            fieldmap_ft = convert(Array{ComplexF32, 5},fieldmap_ft)

            # Saving reshaped calibration scan for reference
            # fieldmap = fftshift(ifft(ifftshift(fieldmap_ft,[1,2]),[1,2]),[1,2])
            # fieldmap = reverse(fieldmap,dims = 1)


            # matwrite(string(params_sv[:path],"acq/fm_",params_sv[:scan],".mat"), Dict("fieldmap" => fieldmap))

            if params_sv[:fmri] == 1
                fm = matread(string("./",params_sv[:directory],"acq/fm_", params_sv[:scan],".mat"))
            else
                fm = matread(string("./",params_sv[:directory],"acq/fm_",params_sv[:scan][1:2],"_01.mat"))
            end
            fm = fm["fieldmap"]
            cs = convert(Array{ComplexF32, 4},cs)
        else
            @info ("Calculating Sensitivity maps ...")
            ### Get Sensitivity maps:
            # fieldmap_ft = matread(string(params_sv[:path],"raw/fieldmap_9ech_b0_",params_sv[:scan][1],params_sv[:scan][4:5],"_1_6_1_6_1_F",".mat"))

            # if params_sv[:fmri] == 1
            #     # fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:scan][1],params_sv[:scan][end-1:end],"_fieldmap.mat"))
            #     fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:scan][1],params_sv[:scan][end-1:end],"_fieldmap.mat"))
            # else
            #     fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:scan][1],"01_fieldmap.mat"))
            # end

            # # Loading fieldmap and saving in image domain
            # fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:fieldmap],"_fieldmap.mat"))

            # # fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:scan],"_fieldmap.mat"))
            # fieldmap_ft = fieldmap_ft["b0"]
            # fieldmap_ft = convert(Array{ComplexF32, 5},fieldmap_ft)

            # # Saving reshaped calibration scan for reference
            # fieldmap = fftshift(ifft(ifftshift(fieldmap_ft,[1,2]),[1,2]),[1,2])
            # fieldmap = reverse(fieldmap,dims = 1)

            # matwrite(string(params_sv[:path],"acq/fm_",params_sv[:scan],".mat"), Dict("fieldmap" => fieldmap))

            lowres_gre = imresize(fieldmap,(params_sv[:nx],params_sv[:ny],params_sv[:sl],size(params_sv[:b0_te])[1]))
            lowres_gre = dropdims(sqrt.(sum(abs.(lowres_gre).^2; dims=5));dims=5)

            # Create Espirit sensitivity maps with MRIReco
            calib_size = 24
            calib_size = calib_size/2
            calibration = fftshift(fft(ifftshift(fieldmap_ft,3),3),3)
            calibration = calibration[Int(end/2+1-calib_size):Int(end/2+calib_size),Int(end/2+1-calib_size):Int(end/2+calib_size),Int(end/2+1-6):Int(end/2+6),1,:]

            cs = espirit(calibration,(params_sv[:nx],params_sv[:ny],params_sv[:sl_no_ov]))

            cs = dropdims(reverse(cs,dims = 1); dims = 5)  # Original

            # cs = dropdims(cs; dims=5)
            # cs = dropdims(reverse(cs,dims = 1); dims = 5)
            # cs = reverse(cs, dims = 3)
            # cs = cs.*-1

            # For now I save the Sensitivity maps and fm in MAT format
            matwrite(string(params_sv[:path],"acq/cs_",file_name,".mat"), Dict("coil_sens" => cs))
            # matwrite(string(params_sv[:path],"acq/fm_",file_name,".mat"), Dict("fieldmap" => fieldmap))
           
        end
    end

    ## Load sensitivity maps mat format (created from BART)
    # cs = matread(string(params_sv[:path],"acq/cs_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],".mat")); cs = cs["coil_sens"];
    # cs = convert(Array{ComplexF32,4},cs);
    # cs = reverse(cs,dims =1 );			# AMM: it looks like I need to reverse the dim, not sure if always

    if params_sv[:do_b0_corr]
    ### Obtain B0 with Fessler approach
    # cs = matread(string(params_sv[:path],"acq/cs_",params_sv[:scan],"_","120","_","122","_",params_sv[:sl],".mat")); cs = cs["coil_sens"];
    # cs = convert(Array{ComplexF32,4},cs);
    # cs ./= sqrt.(sum(abs2, cs; dims=4)) # normalize by SSoS
    
    # mask = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/mask_copy.nii"))
    # mask = mask.raw
    # mask = isone.(mask)
    # fm_mag = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/fieldmap_mag.nii"))
    # fm_mag = fm_mag.raw
    # fm_ph = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/fieldmap_ph.nii"))
    # fm_ph = fm_ph.raw
    # fm = fm_mag.*exp.(1im.*fm_ph)
    # fm = matread(string(params_sv[:path],"raw/",params_sv[:scan],"_fm.mat")); fm = fm["b0"];
    if isfile(string(params_sv[:path],"acq/b0_",file_name,".mat"))
        @info ("Off-resonance map exists... Re-calculate? y/n")
        input = readline()
    else
        input = ""
    end
    if input == "n"
            @info ("Loading Off-resonance map ...")
            b0 = matread(string(params_sv[:path],"acq/b0_",file_name,".mat")); b0 = b0["b0"]
            # b0 = matread("/usr/share/sosp_vaso/data/03032023_sv/tmp/b0_masked.mat"); b0 = b0["b0"]
    else
        @info ("Calculating Off-resonance map ...")
        fieldmap = matread(string(params_sv[:path],"acq/fm_",params_sv[:scan],".mat")); fieldmap = fieldmap["fieldmap"];
        smap = imresize(cs,(size(fieldmap)[1],size(fieldmap)[2],size(fieldmap)[3],size(fieldmap)[5]))

        # smap ./= sqrt.(sum(abs2, smap; dims=4)) # normalize by SSoS
        mask = cs[:,:,:,1]
        mask[abs.(mask) .> 0] .= 1
        mask = imresize(mask,(size(fieldmap)[1],size(fieldmap)[2],size(fieldmap)[3]))
        mask = isone.(mask)
        fm = fieldmap
        fm = permutedims(fm,(1,2,3,5,4))
        fm = fm.*1e3   # Scale data

        # Bias correct..
        params_sv[:b0_te] =  params_sv[:b0_te] * 1e-3 * 1s # echo times in sec
        ######### Temp:
        ydata = fm[:,:,:,:,1:3]; echotime = params_sv[:b0_te][1:3]
        # ydata = fm; echotime = params_sv[:b0_te]

        # smap = cs;
        yik_sos = sum(conj(smap) .* ydata; dims=4) # coil combine
        yik_sos = yik_sos[:,:,:,1,:] # (dims..., ne)

        (yik_sos_scaled, scale) = b0scale(yik_sos, echotime) # todo
        # finit = b0init(ydata, echotime; smap)
        yik_scale = ydata / scale # original
        # yik_scale = ydata

        ## Get the mask from sosp_vaso
        # tmp = abs.(yik_sos[:,:,:,1])
        # mask = tmp.>=maximum(tmp).*0.1

        @infiltrate

        finit = b0init(ydata, echotime; smap)

        fmap_run = (niter, precon, track; kwargs...) ->
            b0map(yik_scale, echotime; finit, smap, mask,
            order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)
        function runner(niter, precon; kwargs...)
            (fmap, times, out) = fmap_run(niter, precon, true; kwargs...)
            # (fmap, _, out) = fmap_run(niter, precon, true; kwargs...) # tracking run
            # (_, times, _) = fmap_run(niter, precon, false; kwargs...) # timing run
            return (fmap, out.fhats, out.costs, times)
        end;
        if !@isdefined(fmap_cg_d)
            niter_cg_d = 400  # 300, 120
            (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag)
            # (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :ichol)
        end

        b0 = fmap_cg_d

        @infiltrate
        
        if params_sv[:field_strength] == 7
            b0 = b0.*2π.*-1    # Original 
        elseif params_sv[:field_strength] == 9
            b0 = b0.*2π.*-1
        elseif params_sv[:field_strength] == 7im
            b0 = b0.*2π.*-1
        end

        @infiltrate
        # b0 = abs.(b0).*-1
        # b0 = reverse(b0, dims=2)
        # b0 = b0.*(1.5)
        
        b0 = ustrip.(b0)
        b0 = imresize(b0,(params_sv[:nx],params_sv[:ny],params_sv[:sl]))
        b0 = 1im.*b0
        b0 = convert(Array{ComplexF32,3},b0)

        @infiltrate

        # Saving NIFTI of GRE first echo
        gre_echo1 = abs.(yik_sos[:,:,:,1])
        gre_echo1 = imresize(gre_echo1,(params_sv[:nx],params_sv[:ny],params_sv[:sl]))
        # gre_echo1 = reverse(gre_echo1,dims = 1)
        # fm = imresize(fm,(params_sv[:nx],params_sv[:ny],params_sv[:sl]))
        gre_echo1 = NIVolume(gre_echo1)
        niwrite(string(params_sv[:path],"tmp/",params_sv[:scan],"_gre_1echo.nii"),gre_echo1)
        # b0[findall(>(500),abs.(b0))] .= Complex(0)
        # b0[findall(<(-500),abs.(b0))] .= Complex(0)

        fn_save_nii(imag.(b0),"B0_jeff")

        mask = imresize(mask,(params_sv[:nx],params_sv[:ny],params_sv[:sl]))

        @infiltrate
        # Masking t2str map
        if params_sv[:do_t2s_corr]
            t2s = matread(string(params_sv[:path],"acq/t2s_",file_name,".mat"));  t2s = t2s["t2str_map"]
            t2s = t2s.*mask
            matwrite(string(params_sv[:path],"acq/t2s_",file_name,".mat"), Dict("t2str_map" => t2s))
        end

        # Saving B0 map and mask, for now only in MAT format
        matwrite(string(params_sv[:path],"acq/b0_",file_name,".mat"), Dict("b0" => b0))
        # matwrite(string(params_sv[:path],"acq/mask_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],".mat"), Dict("mask" => mask))
    
        end
    end

    if params_sv[:do_t2s_corr]
        t2s = matread(string(params_sv[:path],"acq/t2s_",file_name,".mat"));  t2s = t2s["t2str_map"]
        b0 = b0+t2s #(t2s.*1000)
        b0 = convert(Array{ComplexF32,3},b0)
    end

    # reco parameters directory
    params = Dict{Symbol, Any}()
    if params_sv[:do_pi_recon]
        params[:reco] = "multiCoil"; # 'multiCoil'
        params[:senseMaps] = cs
    end
    params[:regularization] = "L1";
    params[:λ] = 1.e-2;  # 1.e-2
    params[:iterations] = 20; # (10)
    params[:solver] = "admm"; # cgnr (L2), fista (L1), admm(L1)
    params[:method] = "nfft"; # nfft, leastsquare

    # set reconstruction size
    if params_sv[:is2d]
        params[:reconSize] = (params_sv[:nx],params_sv[:ny]);
        sli_par = params_sv[:sl];  # original
        # sli_par = 10
    else
        params[:reconSize] = (params_sv[:nx],params_sv[:ny],params_sv[:sl]);
        sli_par = 1;
    end

    contrasts = params_sv[:contrasts];
    # if only want to recon 1 rep, use params[:rep_recon]
    # if params_sv[:multiRepetitions] == false
    #     params_sv[:repetitions] = params_sv[:rep_recon];
    #     f_rep = params_sv[:rep_recon];
    # else
    #     f_rep = 1
    # end

    # Get calibration matrix for B0 dynamic correction
    if params_sv[:drift] == "_drift" && params_sv[:do_b0_corr] == true
        A, B, sh_basis, ΔB0, b_A = getCalibrationMatrix(params_sv,b0,fm)
    end

    # Motion correction, calibration matix
    if params_sv[:mcorr] == "_mCorr"
        @time calib = motionCorrection(params_sv,fm,cs)
    end

    # AMM: Temp, trying to scale B0...
    if params_sv[:do_b0_corr]
        @infiltrate
        b0 = Array{ComplexF32}(b0.*(1.5))
        # b0 = Array{ComplexF32}(b0.*(pi/2))
        # b0 = Array{ComplexF32}(b0.-((20*2*pi)*im))
    end

    j=1:length(contrasts)
    for j=1:length(contrasts)
        # for i=f_rep:params_sv[:repetitions]
        for i = params_sv[:rep_recon][1]:params_sv[:rep_recon][end]
            for k=1:sli_par
                
                if params_sv[:is2d]
                    # Subseting cs and b0    
                    params[:senseMaps] = reshape(cs[:,:,k,:],(params_sv[:nx],params_sv[:ny],1,params_sv[:nCoils]))
                    if params_sv[:do_b0_corr]
                        params[:correctionMap] = reshape(b0[:,:,k],(params_sv[:nx],params_sv[:ny],1))
                    end
                end
                if params_sv[:seq] == 2 || params_sv[:seq] == 4 
                    if params_sv[:id] == "2d"
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_r",i,"_sl",params_sv[:sl_reco],"_",params_sv[:id],"_",params_sv[:traj_type],params_sv[:pdork],params_sv[:rdork],params_sv[:idork] ,".h5"));
                    else
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,params_sv[:pdork],params_sv[:rdork],params_sv[:idork] ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    end
                else
                    if params_sv[:id] == "2d"
                        file =ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_r",i,"_sl",k,"_",params_sv[:id],"_",params_sv[:traj_type],params_sv[:pdork],params_sv[:rdork],params_sv[:idork]  ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    else
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,params_sv[:pdork],params_sv[:rdork],params_sv[:idork]  ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",contrasts[j],"_r",i,"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    end
                end

                @infiltrate
                
                acqData = AcquisitionData(file)

                # # get the number of samples per readout, if its 3d, need to divide by # of slices
                # if params_sv[:is2d]
                #     n_samples = acqData.traj[1].numSamplingPerProfile;
                # else
                #     n_samples = acqData.traj[1].numSamplingPerProfile/params_sv[:sl];
                # end
                n_samples = acqData.traj[1].numSamplingPerProfile;

                # # adjust TE and TAQ after read-in, create times vector for one readout
                # tAQ = (n_samples-1) * params_sv[:dt];
                # acqData.traj[1].AQ = tAQ;                   # Lars: important for B0 correction
                # acqData.traj[1].TE = params_sv[:TE];
                # times = params_sv[:TE] .+ collect(0:params_sv[:dt]:tAQ);
                # acqData.traj[1].times = vec(times).+0.001;

                times = params_sv[:times][:]

                times = times .- params_sv[:TE]

                # if its 3D, repeat the times vector for each slice
                if !params_sv[:is2d]
                    times = repeat(times,Int(floor(params_sv[:sl]/params_sv[:kz])))';
                end

                if params_sv[:mcorr] == "_mCorr" || params_sv[:drift] == "_drift"
                    # Getting the navigator data
                    if !params_sv[:is2d]
                        nav_range = Int(params_sv[:ro_samples]*(params_sv[:sl]/2+1)+1):Int(params_sv[:ro_samples]*(params_sv[:sl]/2+1)+1+300-1)
                    else
                        nav_range = 1:299
                    end
                    nav = (acqData.kdata[1][nav_range,:])
                    nav_times = params_sv[:times][1:299]
                end

                ####### Motion correction
                if params_sv[:mcorr] == "_mCorr"
                    @info ("Motion Correction...")
                    @infiltrate
                    nav_range = findall(nav_times -> nav_times>params_sv[:b0_te][1].*1e-3, nav_times)[1]-20:findall(nav_times -> nav_times>params_sv[:b0_te][1].*1e-3, nav_times)[1]+20
                    nav_mCorr = mean(nav[nav_range,:], dims=1)'
                    nav_mCorr = abs.(nav_mCorr)
                    x_motion = nav_mCorr \ calib
                end

                ####### B0 dynamic off-resonance correction
                if params_sv[:do_b0_corr]
                    b0_drift = Array{ComplexF32}(undef,size(b0))
                    if params_sv[:directory][10:13]=="2022"
                        b0_drift = b0#.*(-1)
                    else
                        b0_drift = b0
                    end

                    if params_sv[:drift] == "_drift"
                    
                        # #### Zero-th order DORK correction.... (not working properly)
                        # # Here I am trying to add the exta Hz difference btw partitions
                        del_omg_v = matread(string(params_sv[:path],"tmp/",params_sv[:scan],"_del_omg_n_v.mat"))
                        del_omg_b = matread(string(params_sv[:path],"tmp/",params_sv[:scan],"_del_omg_n_b.mat"))
                        
                        if params_sv[:contrasts][j] == "v"
                            del_omg = del_omg_v["del_omg_n_v"]
                        elseif params_sv[:contrasts][j] == "b"
                            del_omg = del_omg_b["del_omg_n_b"]
                        end

                        # b0_drift = Array{ComplexF32}(b0_drift.+(2*pi*del_omg[i].*1im))

                        # ##### AMM: Temp, trying to get 
                        # file1=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type]  ,".h5"));
                        # acqData1 = AcquisitionData(file1)
                        # #####

                        #### Higher order correction.... (Wallace approach)
                        # Taking 10 samples centered at TE of B0 (First echo)

                        @infiltrate

                        
                        # tmp = ((params_sv[:b0_te][1].*1e-3)+nav_times[1]) .<= nav_times
                        # tmp = findall(isone,tmp)[1]
                        
                        
                        # nav_range = findall(nav_times -> nav_times>params_sv[:b0_te][tmp-1].*1e-3, nav_times)[tmp-1]-5:findall(nav_times -> nav_times>params_sv[:b0_te][tmp-1].*1e-3, nav_times)[tmp-1]+5
                        nav_range = findall(nav_times -> nav_times>params_sv[:b0_te][1].*1e-3, nav_times)[1]-5:findall(nav_times -> nav_times>params_sv[:b0_te][1].*1e-3, nav_times)[1]+5

                        nav_drift = mean(nav[nav_range,:], dims=1)'

                        # #### Trying to solve the sys of equations with i and r together
                        ###
                        YY = vcat(real.(nav_drift),imag.(nav_drift))
                        AA = vcat(real.(A),imag.(A))
                        b = AA \ YY
                        ####

                        # b = imag.(A) \ (imag.(nav))
                        # b = real.(A) \ (real.(nav))

                        # δB0 = sh_basis * (b)
                        # B is the SH decomposition of the B0 map

                        δB0 = B * (b)
                        δB0 = reshape(δB0,params_sv[:nx],params_sv[:ny],params_sv[:sl])

                        # # AMM: Adding extra Hz to test recon
                        # @infiltrate
                        # δB0 = δB0 .+ (7)

                        mask = deepcopy(b0)
                        mask[abs.(mask) .> 0] .= 1
                        
                        # Extra Hz from Global drift, values from rep DORKs correction
                        δB0 = (δB0.*2*pi).-(abs(del_omg[i]))
                        @infiltrate

                        # Adding del_omg to positive B0 values, and resting from negative...
                        # idx_p = findall(imag.(b0_drift).>0)
                        # idx_n = findall(imag.(b0_drift).<0)
                        # b0_drift[idx_p] = b0_drift[idx_p] .+ ((abs(del_omg[i])*im))
                        # b0_drift[idx_n] = b0_drift[idx_n] .- ((abs(del_omg[i])*im))
                        # @infiltrate
                        # tmp = matread("/usr/share/sosp_vaso/data/05192023_sv/tmp/b0_drift_test.mat")
                        # tmp = matread("/usr/share/sosp_vaso/data/05192023_sv/tmp/b0_drift_test_omg_over_t.mat")
                        # b0_drift = tmp["xx_ft"].*1*im

                        δB0 = δB0.*mask

                        # jim((δB0), "δB0"; color=:jet)

                        # It looks like for VB17 I need to change the sign of the B0map
                        if params_sv[:directory][10:13]=="2022"
                            b0_drift = Array{ComplexF32}(b0_drift-((δB0).*im))
                        else
                            # b0_drift = Array{ComplexF32}(b0_drift+(((δB0.-del_omg[i]).*2*pi).*im))
                            b0_drift = Array{ComplexF32}((b0_drift)-(δB0.*im))
                        end
                        
                        @infiltrate

                    end

                    # AMM: Temp, trying to change b0
                    @infiltrate
                    # b0_drift = Array{ComplexF32}(b0_drift.-((20*2*pi)*im))
                    # rmap = matread(string(params_sv[:path],"tmp/rmap.mat")); rmap = rmap["rmap"];
                    # b0_drift = b0+rmap
                    params[:correctionMap] = b0_drift
                end

                acqData.traj[1].times = vec(times);
                acqData.traj[1].TE = params_sv[:TE];
                acqData.traj[1].AQ = maximum(times);

                @infiltrate

                ######## Reconstruction ###################
                @info ("Starting reconstruction....")
                params = merge(defaultRecoParams(), params);

                @info("Stop... Before recon...")
                @infiltrate

                @time Ireco = reconstruction(acqData, params);
                
                # Removing extra slices due to phase oversampling
                if params_sv[:sl] > params_sv[:sl_no_ov]
                    tmp = (params_sv[:sl]-params_sv[:sl_no_ov])/2
                    Ireco = Ireco[:,:,Int(tmp+1):Int(end-tmp),:,:,:]
                end

                if params_sv[:plt]
                    imshow(abs.(Ireco[:,:,:,1,1]));
                end

                if params_sv[:do_pi_recon]
                    Ireco_mag = NIVolume(abs.(Ireco[:,:,:,1,1,:]))
                else
                    Ireco_mag = sqrt.(sum((abs.(Ireco).^2),dims=5))
                    Ireco_mag = NIVolume(abs.(Ireco_mag[:,:,:,1,:,:]))
                end

                if params_sv[:is2d]
                    save_name = string(params_sv[:path],"recon/2d/",params_sv[:scan],"_",contrasts[j],"_sl",k,"_rep",i,"_",params_sv[:id],"_",params_sv[:traj_type],params_sv[:pdork],params_sv[:rdork],params_sv[:idork],params_sv[:drift] );
                    # save_name = string(params_sv[:path],"recon/2d/",params_sv[:scan],"_",contrasts[j],"_rep",i,"_",params_sv[:id],"_",params_sv[:traj_type] );
                    # save_name = string(params_sv[:path],"recon/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] );
                else
                    save_name = string(params_sv[:path],"recon/3d/",params_sv[:scan],"_",contrasts[j],"_rep",i,"_",params_sv[:id],"_",params_sv[:traj_type],params_sv[:pdork],params_sv[:rdork],params_sv[:idork],params_sv[:drift]);
                    # save_name = string(params_sv[:path],"recon/3d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] );
                end

                if params_sv[:do_b0_corr] 
                    save_name = string(save_name,"_b0_",params_sv[:b0_type])
                end

                if params_sv[:do_t2s_corr] 
                    save_name = string(save_name,"_t2s.nii")
                else
                    save_name = string(save_name,".nii")
                end

                niwrite(save_name,Ireco_mag);

                if params_sv[:save_ph] == 1
                    save_name = string(save_name[1:end-4],"_ph",".nii")
                    Ireco_ph = NIVolume(angle.(Ireco[:,:,:,1,1,:]))
                    niwrite(save_name,Ireco_ph)
                end
                if params_sv[:is2d]
                    @info string("Done reconstructing scan ",params_sv[:scan],  " repetition ",i, " contrast ", contrasts[j], " slice ",k )
                else
                    @info string("done reconstructing scan ",params_sv[:scan], " repetition ", i, " contrast ",  contrasts[j])
                end

            end
        # Trying to clear some memory
        ccall(:malloc_trim, Cvoid, (Cint,), 0) 
        GC.gc()
        end
    end

end

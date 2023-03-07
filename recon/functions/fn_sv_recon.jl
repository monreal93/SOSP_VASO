function fn_sv_recon(params_sv::Dict{Symbol,Any})
    if params_sv[:fmri] == 1
        file_name = string(params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl]);
    else
        file_name = string(params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl]);
    end
    
    if isfile(string(params_sv[:path],"acq/cs_",file_name,".mat"))
        @info ("Sensitivity map exists... Re-calculate? y/n")
        input = readline()
    else
        input = ""
    end
    if input == "n"
        @info ("Loading Sensitivity maps ...")
        cs = matread(string(params_sv[:path],"acq/cs_",file_name,".mat")); cs = cs["coil_sens"];
    else
        @info ("Calculating Sensitivity maps ...")
        ### Get Sensitivity maps:
        # fieldmap_ft = matread(string(params_sv[:path],"raw/fieldmap_9ech_b0_",params_sv[:scan][1],params_sv[:scan][4:5],"_1_6_1_6_1_F",".mat"))

        fieldmap_ft = matread(string(params_sv[:path],"raw/b0_s02_fieldmap.mat"))
        # fieldmap_ft = matread(string(params_sv[:path],"raw/b0_",params_sv[:scan],"_fieldmap.mat"))
        fieldmap_ft = fieldmap_ft["b0"]
        fieldmap_ft = convert(Array{ComplexF32, 5},fieldmap_ft)

        # Saving reshaped calibration scan for reference
        fieldmap = fftshift(ifft(ifftshift(fieldmap_ft,[1,2]),[1,2]),[1,2])
        fieldmap = reverse(fieldmap,dims = 1)
        lowres_gre = imresize(fieldmap,(params_sv[:nx],params_sv[:ny],params_sv[:sl],size(params_sv[:b0_te])[1]))
        lowres_gre = dropdims(sqrt.(sum(abs.(lowres_gre).^2; dims=5));dims=5)

        # Create Espirit sensitivity maps with MRIReco
        calibration = fftshift(fft(ifftshift(fieldmap_ft,3),3),3)
        calibration = calibration[Int(end/2+1-12):Int(end/2+12),Int(end/2+1-12):Int(end/2+12),:,1,:]
        cs = espirit(calibration,(params_sv[:nx],params_sv[:ny],params_sv[:sl]))
        cs = dropdims(reverse(cs,dims = 1); dims = 5)
        # cs = cs.*-1

        # For now I save the Sensitivity maps and fm in MAT format
        matwrite(string(params_sv[:path],"acq/cs_",file_name,".mat"), Dict("coil_sens" => cs))
        matwrite(string(params_sv[:path],"acq/fm_",file_name,".mat"), Dict("fieldmap" => fieldmap))
        # matwrite(string(params_sv[:path],"acq/fm_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],".mat"), Dict("fieldmap" => fieldmap))
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
    else
        @info ("Calculating Off-resonance map ...")
        fieldmap = matread(string(params_sv[:path],"acq/fm_",params_sv[:scan],".mat")); fieldmap = fieldmap["fieldmap"];
        smap = imresize(cs,(size(fieldmap)[1],size(fieldmap)[2],size(fieldmap)[3],size(fieldmap)[5]))
        @infiltrate
        # smap ./= sqrt.(sum(abs2, smap; dims=4)) # normalize by SSoS
        mask = cs[:,:,:,1]
        mask[abs.(mask) .> 0] .= 1
        mask = imresize(mask,(size(fieldmap)[1],size(fieldmap)[2],size(fieldmap)[3]))
        mask = isone.(mask)
        fm = fieldmap
        fm = permutedims(fm,(1,2,3,5,4))
        fm = fm.*1e3   # Scale data
        params_sv[:b0_te] =  params_sv[:b0_te] * 1f-3 * 1s # echo times in sec
        ######### Temp:
        ydata = fm; echotime = params_sv[:b0_te]
        # smap = cs;
        yik_sos = sum(conj(smap) .* ydata; dims=4) # coil combine
        yik_sos = yik_sos[:,:,:,1,:] # (dims..., ne)
        (yik_sos_scaled, scale) = b0scale(yik_sos, echotime) # todo
        # finit = b0init(ydata, echotime; smap)
        @infiltrate
        yik_scale = ydata / scale
        fmap_run = (niter, precon, track; kwargs...) ->
            b0map(yik_scale, echotime; smap, mask,
            order=1, l2b=-4, gamma_type=:PR, niter, precon, track, kwargs...)
        function runner(niter, precon; kwargs...)
            (fmap, times, out) = fmap_run(niter, precon, true; kwargs...)
            # (fmap, _, out) = fmap_run(niter, precon, true; kwargs...) # tracking run
            # (_, times, _) = fmap_run(niter, precon, false; kwargs...) # timing run
            return (fmap, out.fhats, out.costs, times)
        end;
        if !@isdefined(fmap_cg_d)
            niter_cg_d = 300
            (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag)
        end
        fn_save_nii(ustrip.(fmap_cg_d),"B0_jeff")
        b0 = fmap_cg_d
        b0 = b0.*2π.*-1
        b0 = reverse(b0,dims = 1)
        b0 = ustrip.(b0)
        b0 = imresize(b0,(params_sv[:nx],params_sv[:ny],params_sv[:sl]))
        b0 = 1im.*b0;
        b0 = convert(Array{ComplexF32,3},b0)
        mask = imresize(mask,(params_sv[:nx],params_sv[:ny],params_sv[:sl]))
        # Saving B0 map and mask, for now only in MAT format
        matwrite(string(params_sv[:path],"acq/b0_",file_name,".mat"), Dict("b0" => b0))
        # matwrite(string(params_sv[:path],"acq/mask_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],".mat"), Dict("mask" => mask))
        @infiltrate    
        end
    end

    # reco parameters directory
    params = Dict{Symbol, Any}()
    params[:reco] = "multiCoil";
    params[:regularization] = "L1";
    params[:λ] = 1.e-2;  # 1.e-2
    params[:iterations] = 10;
    params[:solver] = "admm"; # cgnr (L2), fista (L1), admm(L1)
    params[:senseMaps] = cs;
    params[:method] = "nfft"; # nfft, leastsquare

    # set reconstruction size
    if params_sv[:is2d]
        params[:reconSize] = (params_sv[:nx],params_sv[:ny]);
        sli_par = params_sv[:sl];  # original
        sli_par = 10
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
                if params_sv[:seq] == 2
                    if params_sv[:id] == "2d"
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_r",i,"_sl",params_sv[:sl_reco],"_",params_sv[:id],"_",params_sv[:traj_type],params_sv[:dork] ,".h5"));
                    else
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,params_sv[:dork] ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    end
                else
                    if params_sv[:id] == "2d"
                        file =ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_r",i,"_sl",k,"_",params_sv[:id],"_",params_sv[:traj_type],params_sv[:dork]  ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    else
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,params_sv[:dork]  ,".h5"));
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
                    times = repeat(times,Int(params_sv[:sl]/params_sv[:kz]))';
                end

                # Add B0 map and correct for scanner drift
                if params_sv[:do_b0_corr]
                    b0_drift = Array{ComplexF64}(undef,size(b0))
                    b0_drift = b0
                    # b0_drift[findall(!iszero,b0_drift)] = b0_drift[findall(!iszero,b0_drift)].-(1im.*0.14*f_rep*2*pi)
                    # b0_drift[findall(!iszero,b0_drift)] = b0_drift[findall(!iszero,b0_drift)].-(1im.*0.14*72*2*pi)
                    params[:correctionMap] = b0_drift
                end

                acqData.traj[1].times = vec(times);
                acqData.traj[1].TE = params_sv[:TE];
                acqData.traj[1].AQ = maximum(times);

                ######## Reconstruction ###################
                @info ("Starting reconstruction....")
                params = merge(defaultRecoParams(), params);
                @time Ireco = reconstruction(acqData, params);

                if params_sv[:plt]
                    imshow(abs.(Ireco[:,:,:,1,1]));
                end

                Ireco_mag = NIVolume(abs.(Ireco[:,:,:,1,1,:]));


                if params_sv[:is2d]
                    save_name = string(params_sv[:path],"recon/2d/",params_sv[:scan],"_",contrasts[j],"_sl",k,"_rep",i,"_",params_sv[:id],"_",params_sv[:traj_type] );
                    # save_name = string(params_sv[:path],"recon/2d/",params_sv[:scan],"_",contrasts[j],"_rep",i,"_",params_sv[:id],"_",params_sv[:traj_type] );
                    # save_name = string(params_sv[:path],"recon/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] );
                else
                    save_name = string(params_sv[:path],"recon/3d/",params_sv[:scan],"_",contrasts[j],"_rep",i,"_",params_sv[:id],"_",params_sv[:traj_type] );
                    # save_name = string(params_sv[:path],"recon/3d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] );
                end

                if params_sv[:do_b0_corr] 
                    save_name = string(save_name,"_b0_",params_sv[:b0_type],"_mrreco.nii")
                else
                    save_name = string(save_name,"_mrreco.nii")
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
        end
    end

end

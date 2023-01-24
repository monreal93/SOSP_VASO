function fn_sv_recon(params_sv::Dict{Symbol,Any})
    # load sensitivity maps mat format
    cs = matread(string(params_sv[:path],"acq/cs_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],".mat")); cs = cs["coil_sens"];

    # converting to ComplexF32, that is what MRIReco.jl expects
    cs = convert(Array{ComplexF32,4},cs);
    cs = reverse(cs,dims =1 );			# AMM: it looks like I need to reverse the dim, not sure if always

    # load B0 map, romeo/skope-i
    if params_sv[:do_b0_corr]   
        if params_sv[:b0_type] == "romeo"
            b0 = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/B0_sm.nii")) # From Romeo
            b0 = b0.raw;
            replace!(b0, NaN=>0);
            # AMM: Temp: Trying to scale
            b0 = b0.*2π.*-1
            # b0 = b0.*-1 # (-2.5), it seems like I have to adjust this value "scale" each b0 differently...
            b0 = reverse(b0,dims = 1);
        elseif params_sv[:b0_type] == "gilad"
            b0 = niread(string(params_sv[:path],"acq/gilad/b0_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"_gilad.nii")); # From Gilad's
            b0 = b0.raw;
            # AMM: Temp: Trying to scale
            b0 = b0.*π;
            b0 = b0.*-1;
            b0 = reverse(b0,dims = 2);
        elseif params_sv[:b0_type] == "skope"
            b0 = niread(string(params_sv[:path],"acq/skope/b0_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"_skope.nii")); # From Skope-i
            b0 = b0.raw;
            # AMM: Temp: Trying to scale
            # b0 = b0.*π;
            # b0 = b0.*-1;
            b0 = b0.*-0.75;
            b0 = reverse(b0,dims = 1);
        elseif params_sv[:b0_type] == "fessler"
            b0 = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/B0_jeff_80iter_o1.nii")); 
            b0 = b0.raw;
            replace!(b0, NaN=>0);
            # AMM: Temp: Trying to scale
            b0 = b0.*2π.*-1
            # b0 = b0.*-1 # (-2.5), it seems like I have to adjust this value "scale" each b0 differently...
            b0 = reverse(b0,dims = 1);
        end

        b0 = 1im.*b0;
        b0 = convert(Array{ComplexF32,3},b0);
    end

    # ### Obtain B0 with Fessler approach
    # cs = matread(string(params_sv[:path],"acq/cs_",params_sv[:scan],"_","120","_","122","_",params_sv[:sl],".mat")); cs = cs["coil_sens"];
    # cs = convert(Array{ComplexF32,4},cs);
    # # cs ./= sqrt.(sum(abs2, cs; dims=4)) # normalize by SSoS
    # mask = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/mask_copy.nii"))
    # # mask = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/mask.nii"))
    # mask = mask.raw
    # mask = isone.(mask)
    # fm_mag = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/fieldmap_mag.nii"))
    # fm_mag = fm_mag.raw
    # fm_ph = niread(string(params_sv[:path],"acq/romeo/b0_",params_sv[:scan],"_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/fieldmap_ph.nii"))
    # fm_ph = fm_ph.raw
    # fm = fm_mag.*exp.(1im.*fm_ph)
    # # fm = matread(string(params_sv[:path],"raw/",params_sv[:scan],"_fm.mat")); fm = fm["b0"];
    # fm = permutedims(fm,(1,2,3,5,4))
    # fm = fm.*1e3;   # Scale data
    # params_sv[:b0_te] =  params_sv[:b0_te] * 1f-3 * 1s # echo times in sec
    # @infiltrate
    # ######### Temp:
    # ydata = fm; smap = cs; echotime = params_sv[:b0_te]
    # yik_sos = sum(conj(smap) .* ydata; dims=4) # coil combine
    # yik_sos = yik_sos[:,:,:,1,:] # (dims..., ne)
    # (yik_sos_scaled, scale) = b0scale(yik_sos, echotime) # todo
    # finit = b0init(ydata, echotime; smap)
    # yik_scale = ydata / scale
    # fmap_run = (niter, precon, track; kwargs...) ->
    #     b0map(yik_scale, echotime; smap, mask,
    #     order=2, l2b=-4, gamma_type=:PR, niter, precon, track, kwargs...)
    # function runner(niter, precon; kwargs...)
    #     (fmap, _, out) = fmap_run(niter, precon, true; kwargs...) # tracking run
    #     (_, times, _) = fmap_run(niter, precon, false; kwargs...) # timing run
    #     return (fmap, out.fhats, out.costs, times)
    # end;
    # if !@isdefined(fmap_cg_d)
    #     niter_cg_d = 80
    #     (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag)
    # end
    # fn_save_nii(ustrip.(fmap_cg_d),"B0_jeff")
    # @infiltrate
    # ###########

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

    if params_sv[:do_b0_corr]
        params[:correctionMap] = b0;
    end

    contrasts = params_sv[:contrasts];
    # if only want to recon 1 rep, use params[:rep_recon]
    if params_sv[:multiRepetitions] == false
        params_sv[:repetitions] = params_sv[:rep_recon];
        f_rep = params_sv[:rep_recon];
    else
        f_rep = 1
    end

    j=1:length(contrasts)
    for j=1:length(contrasts)
        for i=f_rep:params_sv[:repetitions]
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
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_r",i,"_sl",params_sv[:sl_reco],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    else
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,".h5"));
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    end
                else
                    if params_sv[:id] == "2d"
                    file =ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_r",i,"_sl",k,"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_r",i,"_",params_sv[:traj_type] ,".h5"));
                    # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/2d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    else
                        # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",contrasts[j],"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                        file=ISMRMRDFile(string(params_sv[:path],"ismrmd/3d/",params_sv[:scan],"_",contrasts[j],"_r",i,"_",params_sv[:id],"_",params_sv[:traj_type] ,".h5"));
                    end
                end

                acqData = AcquisitionData(file);

                # get the number of samples per readout, if its 3d, need to divide by # of slices
                if params_sv[:is2d]
                    n_samples = acqData.traj[1].numSamplingPerProfile;
                else
                    n_samples = acqData.traj[1].numSamplingPerProfile/params_sv[:sl];
                end

                # # adjust TE and TAQ after read-in, create times vector for one readout
                # tAQ = (n_samples-1) * params_sv[:dt];
                # acqData.traj[1].AQ = tAQ;                   # Lars: important for B0 correction
                # acqData.traj[1].TE = params_sv[:TE];
                # times = params_sv[:TE] .+ collect(0:params_sv[:dt]:tAQ);
                # acqData.traj[1].times = vec(times).+0.001;

                times = params_sv[:times]

                # if its 3D, repeat the times vector for each slice
                if !params_sv[:is2d]
                    times = repeat(times,Int(params_sv[:sl]/params_sv[:kz]))';
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

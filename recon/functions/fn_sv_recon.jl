###########     To Do
#   1) Check Temp: line 43, do I want to shift it?
#   2) Decide if I want cs in mat file or nifti
#   3)
#   4)
#   5)
###########
function fn_sv_recon(params_sv::Dict{Symbol,Any})
    # file=ISMRMRDFile(string(params_sv[:path],"ismrmd/",params_sv[:scan],"_",params_sv[:id],".h5"));
    # acqData = AcquisitionData(file);

    # Loading data
    ### cs mat format
    cs = matread(string(params_sv[:path],"acq/cs_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],".mat")); cs = cs["coil_sens"];
    ### cs nifti format
    # tmp = niread("/home/amonreal/Documents/PhD/PhD_2021/data/cs_134_134_24_r_msk.nii"); tmp = tmp.raw;
    # itmp = niread("/home/amonreal/Documents/PhD/PhD_2021/data/cs_134_134_24_i_msk.nii"); itmp = itmp.raw;
    # cs = tmp+(itmp*im);
    # Shifting CS slices to make match the gre to the spiral acq
    # cs = circshift(cs, (0,0,-1,0));
    # cs = circshift(cs, (1,0,0,0));

    if params_sv[:is2d]
        cs = cs[:,:,params_sv[:sl_reco],:];
        cs = reshape(cs,(params_sv[:nx],params_sv[:ny],1,params_sv[:nCoils]));
    end
    # converting to ComplexF32, that is what MRIReco.jl expects
    cs = convert(Array{ComplexF32,4},cs);
    cs = reverse(cs,dims =1 );			# AMM: it looks like I need to reverse the dim, not sure if always

    # load B0 map, romeo/skope-i
    if params_sv[:do_b0_corr]   
        if params_sv[:b0_type] == "romeo"
            b0 = niread(string(params_sv[:path],"acq/romeo/",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/B0_masked_sm.nii")); # From Romeo
            # b0 = niread(string(params_sv[:path],"acq/romeo/",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"/b0_app2.nii")); # From Romeo
            b0 = b0.raw;
            replace!(b0, NaN=>0);
            b0 = b0.*2π;
            # AMM: Temp: Trying to scale
            b0 = b0.*1;
            b0 = reverse(b0,dims = 1);
        elseif params_sv[:b0_type] == "gilad"
            b0 = niread(string(params_sv[:path],"acq/gilad/b0_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"_gilad.nii")); # From Gilad's
            b0 = b0.raw;
            b0 = b0.*2π;
            # AMM: Temp: Trying to scale
            # b0 = b0.*π;
            b0 = b0.*0.8;
            b0 = reverse(b0,dims = 2);
        elseif params_sv[:b0_type] == "skope"
            b0 = niread(string(params_sv[:path],"acq/skope-i/b0_",params_sv[:nx],"_",params_sv[:ny],"_",params_sv[:sl],"_skope.nii")); # From Skope-i
            b0 = b0.raw;
            # AMM: Temp: Trying to scale
            b0 = b0.*π;
        end
        # Temp: Do I really want to shift it??? Shifting B0 slices to make match the gre to the spiral acq
        # b0 = circshift(b0, (0,0,1));

        b0 = 1im.*b0;
        b0 = convert(Array{ComplexF32,3},b0);
        # b0 = reverse(b0,dims = 1);			# AMM: it looks like I need to reverse the dim, not sure if always
        # b0 = circshift(b0, (1,0,0));

        if params_sv[:is2d]
            b0 = b0[:,:,params_sv[:sl_reco]];
            b0 = reshape(b0,(params_sv[:nx],params_sv[:ny],1));
        end

        # AMM: Infiltrate
        @infiltrate
    end

    # # get the number of samples per readout, if its 3d, need to divide by # of slices
    # if params_sv[:is2d]
    #     n_samples = acqData.traj[1].numSamplingPerProfile;
    # else
    #     n_samples = acqData.traj[1].numSamplingPerProfile/params_sv[:sl];
    # end

    # # adjust TE and TAQ after read-in, create times vector for one readout
    # tAQ = (n_samples-1) * params_sv[:dt];
    # acqData.traj[1].AQ = tAQ;                   # Lars: important for B0 correction
    # acqData.traj[1].TE = params_sv[:TE];
    # times = params_sv[:TE] .+ collect(0:params_sv[:dt]:tAQ);

    # # if its 3D, repeat the times vector for each slice
    # if !params_sv[:is2d]
    #     times = repeat(times,params_sv[:sl]);
    # end

    # acqData.traj[1].times = times;

    # reco parameters directory
    params = Dict{Symbol, Any}()
    params[:reco] = "multiCoil";

    # set reconstruction size
    if params_sv[:is2d]
        params[:reconSize] = (params_sv[:nx],params_sv[:ny]);
    else
        params[:reconSize] = (params_sv[:nx],params_sv[:ny],params_sv[:sl]);
    end
    params[:regularization] = "L2";
    params[:λ] = 1.e-2;
    params[:iterations] = 10;
    params[:solver] = "cgnr";
    params[:senseMaps] = cs;

    if params_sv[:do_b0_corr]
        params[:correctionMap] = b0;
    end

    # Loop for vaso and bold, also for each echo..
    # contrasts = ["v","b"];
    contrasts = params_sv[:contrasts];
    # if only want to recon 1 rep (it will do rep 2)
    if params_sv[:multiRepetitions] == false
        params_sv[:repetitions] = 2;
        f_rep = 2;
    else
        f_rep = 1;
    end
    # AMM: Infiltrate
    # @infiltrate

    for i=f_rep:params_sv[:repetitions]
        for j=1:length(contrasts)

            if params_sv[:seq] == 2
                file=ISMRMRDFile(string(params_sv[:path],"ismrmd/",params_sv[:scan],"_r",i,"_",params_sv[:id],".h5"));
            else
                file=ISMRMRDFile(string(params_sv[:path],"ismrmd/",params_sv[:scan],"_",contrasts[j],"_r",i,"_",params_sv[:id],".h5"));
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

            times = params_sv[:times];

            # if its 3D, repeat the times vector for each slice
            if !params_sv[:is2d]
                times = repeat(times,Int(params_sv[:sl]/params_sv[:kz]))';
            end

            acqData.traj[1].times = vec(times);
            acqData.traj[1].TE = params_sv[:TE];
            acqData.traj[1].AQ = maximum(times);
            
            ######## Reconstruction ###################
            params = merge(defaultRecoParams(), params);
            @time Ireco = reconstruction(acqData, params);
            Ireco = Ireco.*1e6;

            if params_sv[:plt]
                imshow(abs.(Ireco[:,:,:,1,1]));
            end

            Ireco_nii = NIVolume(abs.(Ireco[:,:,:,1,1]));

            if params_sv[:is2d]
                save_name = string(params_sv[:path],"recon/",params_sv[:scan],"_",contrasts[j],"_sl",params_sv[:sl_reco],"_rep",i,"_",params_sv[:id]);
            else
                save_name = string(params_sv[:path],"recon/",params_sv[:scan],"_",contrasts[j],"_rep",i,"_",params_sv[:id]);
            end

            if params_sv[:do_b0_corr] 
                save_name = string(save_name,"_b0_",params_sv[:b0_type],"_mrreco.nii")
            else
                save_name = string(save_name,"_mrreco.nii")
            end

            niwrite(save_name,Ireco_nii);

            @info string("Done reconstructing ", contrasts[j], " repetition",i)

        end
    end

end

function params = prepare_add_parameters(ks_traj,gx,rf,adc,te0,te1,tr0,tr1,params)

    % Update resolution and mtx size with effective resolution
    rx = 1./(abs(max(ks_traj.kx(:,1)))+abs(min(ks_traj.kx(:,1))));
    if params.gen.ro_type == 'c'
        ry = 1./(abs(max(ks_traj.ky(:,1)))*2);
    else
        ry = 1./(abs(max(ks_traj.ky(:,1)))+abs(min(ks_traj.ky(:,1))));
    end
    % If scan is 2D
    if params.gen.n(3) == 1
        rz = params.gen.res(3);
    else
        rz = 1./(abs(max(ks_traj.kz(:)))+abs(min(ks_traj.kz(:))));
    end
    params.gen.res = [rx,ry,rz];

    params.gen.n(1:2) = round(params.gen.fov(1:2)./params.gen.res(1:2));
    tmp = mod(params.gen.n(1:2),2); 
%     params.gen.n(1:2) = params.gen.n(1:2)+tmp; % making in-plane even
%     params.gen.n(1:2) = max(params.gen.n(1),params.gen.n(2)); % Making it square??
    % Adjust for phase oversampling
    params.gen.n_ov(1:2) = params.gen.n(1:2);
    % ToDo: Check what value I really want in n_ov (to be used in recon)
%     params.gen.n_ov(3) = params.gen.n_ov(3).*params.gen.kz;
    [params.gen.n, params.gen.n_ov] = deal(params.gen.n_ov, params.gen.n);
    
    % converting FA back to degrees
    params.gen.fa = real(params.gen.fa)*180/pi;
    
    %% Adding ro_samples, acqTR,volTR, effTR, TE and other params
    params.gen.ro_samples = adc.numSamples;
    params.gen.TE = te1-te0+(mr.calcDuration(rf)/2);
    params.gen.acqTR = tr1-tr0;
    params.gen.volTR = params.gen.acqTR.*params.gen.n_ov(3).*params.spi.interl;
    % params.gen.effTR = seq_t1-seq_t0;
    params.gen.adc_dwell = adc.dwell;
    if params.gen.ro_type == 's'
        params.gen.ro_time = adc.duration;
    elseif params.gen.ro_type == 'c'
        params.gen.ro_time = adc.duration*params.epi.n_lines;
    end
    % Calculating TI1 for VASO sequences...
    if params.gen.kz_enc == 0
        if params.gen.seq == 1 && params.vaso.bold_ref == 1
            params.vaso.ti1 = ((params.gen.effTR-params.vaso.f_v_delay)/4)+params.vaso.f_v_delay;
        elseif params.gen.seq == 1 && params.vaso.bold_ref == 0
            params.vaso.ti1 = ((params.gen.effTR-params.vaso.f_v_delay)/2)+params.vaso.f_v_delay;
        end
    elseif params.gen.kz_enc == 1
        % ToDo: for a more precise number, add the time for fatsat and
        % ph-enc grad...
        params.vaso.ti1 = params.vaso.f_v_delay;
    end
    if params.spi.type ==3; params.gen.echos = 2; end       % If IN-OUT.. 2 echos

    % Adjusting BW, with correct dwell
    params.spi.bw = 1/params.gen.adc_dwell;
    
    %% Calculating time vector for MRIReco.jl
    % AMM: Here I might need to take into account the gradient delays...
    % julia_time = params.gen.TE+adc.delay+(adc.dwell:adc.dwell:adc.duration);
    if params.gen.ro_type == 's'
        julia_time = repmat(params.gen.TE+(0:adc.dwell:adc.duration-adc.dwell),1,params.spi.interl);
    elseif params.gen.ro_type == 'c'
        julia_time = (0:adc.dwell:(adc.duration*params.epi.n_lines/params.epi.seg)-adc.dwell);
        julia_time = params.gen.TE+ julia_time + (mr.calcDuration(gx)*3) + mr.calcDuration(rf)/2 + 0.6e-3 + 4e-4; % Calibration lines + ~0.6e-3 of gy_pre + 4e-4 of reph-rf_grad
        julia_time = repmat(julia_time,1,params.epi.seg);
    end
    params.gen.t_vector = julia_time;
    
    % Adjusting for different readout types...
    if params.gen.ro_type == 's' && params.spi.type == 3
        params.gen.TE = params.gen.TE+(mr.calcDuration(gx)/2);
        params.gen.ro_samples = params.gen.ro_samples/2;
    end

end

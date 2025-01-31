clear all; clc
% Change directory to path where sosp_vaso git has been cloned
cd /home/amonreal/Documents/PhD/PhD_2025/sosp_vaso/

%% Adding paths 
addpath(genpath("./pulseq/functions"))                                                                  % Functions to create pulseq sequence
addpath(genpath("/home/amonreal/Documents/PhD/tools/pulseq_1.4.1/"))                                    % Pulseq toolbox
addpath(genpath("/home/amonreal/Documents/PhD/tools/tOptGrad_V0.2/minTimeGradient/mex_interface/"))     % Minminue time gradieng Lusitg
addpath(genpath("/home/amonreal/Documents/PhD/tools/pns_prediction/"))                                  % PNS prediction
addpath(genpath("/home/amonreal/Documents/PhD/tools/check_grad_idea_Des/"))                             % Check Forbidden Fq
warning('OFF', 'mr:restoreShape')

%% Define parameters
folder_name = 'simulations_vaso_only';                % Folder will be created in ./data
seq_name = 'sample';                % use sv/abc/sb_n (n for the diff scans at each day)
params.gen.seq = 4;                 % 1-VASO 2-ABC 3-Multi-Echo 4-BOLD
params.gen.field_strength = 9;      % Field Strength (7=7T,7i=7T-impuse_grad,9=9.4T,11=11.7T)
params.gen.pns_check = 1;           % PNS check, .acs files to be added in ./tools/pns_check/grad_files

%%%% General parameters
params.gen.fov = [140 140 24].*1e-3;% FOV (Josh's project 160
params.gen.res = [0.6 0.6 0.6].*1e-3;% Nominal resolution
params.gen.fa = 0;                  % FA in degrees. Set to 0, to use Ernst Angle 
params.gen.vfa = 0;                 % Variable FA (WIP)
params.gen.vfa_cutoff = 33;         % Variable FA cut-off (WIP) in degrees
params.gen.ernst_t1 = 2100e-3;      % T1 to calculate Ernst Angle (s) 7T=(WM-1220e-3)(GM-1800e-3)(blood=2587e-3) 9T=(1425e-3)(2100e-3)
params.gen.te = 6e-3;               % Set to 0 to shortest TE possible
params.gen.multi_te = [12e-3 6e-3]; % Echo times for ME seq=3 start with longest (WIP)
params.gen.tr_delay = 0e-3;         % Delay between acquisitions in sec (0e-3)
params.gen.ro_type = 's';           % 's'-Spiral, 'c'-Cartesiaen (WIP)
params.gen.kz = 1;                  % Acceleration in Kz (1)
params.gen.kz_enc = 0;              % k-space partition encoding 0=linear,1=center-out For Cartesian now only linear encoding
params.gen.pf = 1;                  % Partial fourier in Kz
params.gen.fat_sat = 1;             % Fat saturation (1=yes,0=no)
params.gen.fs_angle = 0;            % Fat sat angle (0=default)
params.gen.fs_interl = 1;           % Fat sat in every shot/interl? (WIP)
params.gen.skope = 0;               % Add skope sync scan and triggers, 0=N0, 1=sep scan, 2=concurrent(center partition), 3=1&2
params.gen.skope_sync = 0;          % Skope pre-scans, only added in skope seq
params.gen.dork = 0;                % extra adc's for DORK correction (WIP)
params.gen.interl_enc = 0;          % Interleave encoding, 0=1-interl 1-plane -> 2-interl 1-plane, 1=1-inter 1-plane -> 1-interl 2-plane
params.gen.ph_oversampling = 0;     % Partition phase oversampling in %, to avoid partition phase-wrap (default=0)
params.gen.echos = 1;               % Echos per RF pulse (WIP)
params.gen.fid_nav = 1;             % FID navigators before each readout for Spiral
params.gen.kz_caipi = 1;            % CAIPI shift in partition direction (1-No,2-half) (WIP)
%%%% ME-GRE parameters
params.gen.me_gre = 0;              % ME GRE calibration scan, 0=NO, 1=same seq, 2=separate seq
params.gen.me_gre_echos = 3;        % ME GRE calibration scan echos (>1)
params.gen.me_gre_tr = 30e-3;       % ME GRE TR (60e-3)
params.gen.me_gre_interl = 42;      % ME GRE Shots (32)

% Spiral parameters
params.spi.type = 0;                % spiral type 0=spiral-Out , 1=spiral-In, 3=In-Out (WIP), 4=In-Out kspace interleavead (WIP)
params.spi.in_out_order = 0;        % 0=In-Out same k-space path (separate vol.), 1=In-Out k-space path shift (WIP)
params.spi.rotate = 'none';         % Spiral rotation ('none','golden','180','120')
params.spi.increment = 'linear';    % Spiral increment mode (for now only 'linear') (WIP)
params.spi.max_grad  = 50;          % Peak gradient amplitude for spiral (mT/m)
params.spi.max_sr = 300;            % Max gradient slew rate for spiral (mT/m/ms) (25for ABC).
params.spi.interl = 2;              % Spiral interleaves
params.spi.vd = 1.3;                % Variability density (accepts negative values)
params.spi.rxy = 3.3;               % In-plane (radial) undersampling
params.spi.rxy_az = 1;              % In-plane (azimuthal) undersampling (WIP)
params.spi.bw = 0e3;                % Spiral BW in Hz 0=Nyquist, highest=1000e3
params.spi.safe_sp = 0;             % Safe spirals (WIP)

% MT pulse parameters
params.mt.mt = 0;                   % Add MT pulse, 0 for reference scan without MT
params.mt.mt_prep = 0;              % MT (pre) dummy scans
params.mt.rf_spoil = 0;             % RF spoliing
params.mt.alpha = 225;              % default=225
params.mt.delta = -650;             % default=-650
params.mt.trf = 0.004;
params.mt.mt_rep = 100e-3;          % How often to play the MT pulse
params.mt.bold = 0;                 % Get BOLD reference acq after MT one (WIP)

% EPI parameters
params.epi.ry = 3;                  % In-plane undersampling
params.epi.pf = 6/8;                % In-plane PF, (1,7/8,6/8)
params.epi.te = [33.6 36.6 38.6 40]*1e-3+0.07; % Echo times for field map
params.epi.seg = 2;                 % EPI Segments
params.epi.tr_delay = 0;            % Delay after each TR, needed to reduce SAR
params.epi.bw_px = 1512;            % BW in Hz (1096)

% VASO parameters
params.vaso.foci = 1;               % FOCI inversion?
params.vaso.bold_ref = 1;           % BOLD reference volume
params.vaso.foci_ampl = 270;        % FOCI amplitude (270)
params.vaso.foci_dur = 10410e-6;    % FOCI duration (10410e-6)
params.vaso.foci_bw = 150;          % FOCI Bandwidth (150)
params.vaso.f_v_delay = 600e-3;     % FOCI-VASO delay (600e-3)
params.vaso.v_b_delay = 10e-3;      % VASO-BOLD delay (10e-3)
params.vaso.b_f_delay = 5e-3;       % BOLD-FOCI delay (5e-3)

%% Set system limits
params = prepare_system_limits(params);

%% Calculations and restrictions
[params,ro_blocks] = prepare_fix_parameters(params);

%% Create blocks
% TR-FOCI
if params.gen.seq == 1
    [B1_foci,phase,rf_complex,Gz_foci,fa] = tr_foci(params);
    rf_foci = mr.makeArbitraryRf(rf_complex,fa*pi/180, 'system', params.gen.lims);%, 'Delay',2e-4);
end

% MT pulse
if params.gen.seq == 2
    MT = mr.makeGaussPulse(params.mt.alpha*pi/180,params.gen.lims,'Duration',params.mt.trf,'FreqOffset',params.mt.delta);
end

% Kz blips
gz_blips = prepare_gz_blips(params);

% Gradient spoilers
grad_area = 1/params.gen.fov(1)*params.gen.n(1)/2;
[gx_spoil,gy_spoil,gz_spoil] = prepare_spoilers(params,grad_area);

% Fat sat
if params.gen.fat_sat || params.gen.me_gre ~= 0
    sat_ppm = -3.45;
    [rf_fs,gx_fs,gy_fs,gz_fs] = prepare_fat_sat(params, sat_ppm);
    grad_area = abs(gx_spoil.area)+abs(gx_fs.area);
    [gx_spoil,gy_spoil,gz_spoil] = prepare_spoilers(params,grad_area);
end


%% Preparing readout elements
[rf_phase_offset,adc_phase_offset] = rf_adc_phase(params,0);
if params.gen.ro_type == 's'
    if params.spi.safe_sp
        [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc_safe_sp(params);
    else
        [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc(params);
    end
    if params.spi.type == 0
        [gx,gy,~,~,adc,params] =  create_spirals_adc_pulseq(spiral_grad_shape,adcSamples,adcDwell,params);
    elseif params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
        [gx,gy,gx_pre,gy_pre,adc,params] =  create_spirals_adc_pulseq(spiral_grad_shape,adcSamples,adcDwell,params);
    end
elseif params.gen.ro_type == 'c'
    [gx,gy_blips,gx_pre,gy_pre,gy_blip_up,gy_blip_down,adc,params]  = create_epi_adc_pulseq(params);
end

% Flip angle, 3.5e-3 (approx rf + rephasing time)
tr_tmp = 3.5e-3+mr.calcDuration(gz_blips(1))+mr.calcDuration(mr.calcDuration(gx)*params.gen.echos)+mr.calcDuration(gx_spoil)+params.gen.tr_delay; % aprox TR
if params.gen.seq == 3
    tr_tmp = tr_tmp +params.gen.multi_te(1);
else
    tr_tmp = tr_tmp +params.gen.te;
end
if and(params.gen.ro_type=='s',or(params.spi.type == 1,params.spi.type == 2)); tr_tmp=tr_tmp+mr.calcDuration(gx_pre(1));end 
params = prepare_flip_angle(tr_tmp,gx, params);

% Create RF and Gz
[rf0,gz,gzReph,params] =  create_rf_gz_pulseq(params,tr_tmp);

% Create FID navigators
if params.gen.fid_nav
    fid_nav_adc = mr.makeAdc(100,params.gen.lims,'Dwell',adc.dwell);
    fid_nav_delay1 = mr.makeDelay(mr.calcDuration(gz_blips(1))-mr.calcDuration(fid_nav_adc)/2); % Delay to center the first FID and ME-GRE TE...
    fid_nav_delay2 = mr.makeDelay(mr.calcDuration(fid_nav_adc)*2);                              % FID delay is twice the FID duration...
end

%% Prepare Fieldmap GRE-ME elements
if params.gen.me_gre > 0 && params.gen.ro_type == 's'
    params_me_gre = params;
    params_me_gre.gen.fat_sat = 1;          % Let's always use fat sat for ME-GRE
%     params_me_gre.gen.res(1:2) = params_me_gre.gen.res(1:2).*1.75; % (2)
    params_me_gre.gen.res(1:2) = 1e-3;  % Fix it at 1.2 mm
    params_me_gre.spi.interl = params.gen.me_gre_interl;  % (32)
    params_me_gre.gen.seg = params_me_gre.spi.interl;
    params_me_gre.spi.rxy = 1;
    params_me_gre.spi.vd = 1;
    params_me_gre.spi.safe_sp = 0;
    params_me_gre.gen.kz = 1;
    params_me_gre.gen.kz_caipi = 1;
    params_me_gre.gen.fa = params_me_gre.gen.fa(1)/(pi/180); % Same as fMRI
    params_me_gre.gen.echos = params.gen.me_gre_echos;
    params_me_gre.spi.rotate = 'none';
%     params_me_gre.gen.ph_oversampling = 0;
    
    %%%%% Original....
    % No under-sampling for GRE-ME:
    params_me_gre.gen.n(3) = round(params_me_gre.gen.fov(3)./params_me_gre.gen.res(3));
    params_me_gre.gen.del_k = 1./(params_me_gre.gen.fov);

%     %%%%% Temp: Increasing FOV a bit... for motion correction simulations
%     params_me_gre.gen.fov(1:2) = params_me_gre.gen.fov(1:2)+10e-3;
%     params_me_gre.gen.fov(3) = params_me_gre.gen.fov(3)+4e-3;
%     params_me_gre.gen.n = round(params_me_gre.gen.fov./params_me_gre.gen.res);
%     tmp = mod(params_me_gre.gen.n(3),2); 
%     params_me_gre.gen.n(3) = params_me_gre.gen.n(3)+tmp; % making in-plane even
%     params_me_gre.gen.n_ov = params_me_gre.gen.n;
%     % No under-sampling for GRE-ME:
%     params_me_gre.gen.del_k = 1./(params_me_gre.gen.fov);

    % GRE-ME spiral MaxG and SR, might need to change if PNS errors
    if params.gen.field_strength == 9
        params_me_gre.spi.max_sr = 120; % (7T-60)
        params_me_gre.spi.max_grad = 40; % (7T-20)
    else
        params_me_gre.spi.max_sr = 60; % (7T-60)
        params_me_gre.spi.max_grad = 20; % (7T-20)
    end
    params_me_gre.spi.rotate = 'none';
    params_me_gre.spi.type = 0;
    [rf_phase_offset_me_gre,adc_phase_offset_me_gre] = rf_adc_phase(params_me_gre,0);
    if params.gen.ro_type == 's'
            [spiral_grad_shape_me_gre,adcSamples_me_gre,adcDwell_me_gre,params_me_gre] = prepare_spirals_rf_grad_adc(params_me_gre);
            [gx_me_gre,gy_me_gre,~,~,adc_me_gre,params_me_gre] =  create_spirals_adc_pulseq(spiral_grad_shape_me_gre,adcSamples_me_gre,adcDwell_me_gre,params_me_gre);
    end
    params_me_gre.gen.ro_samples = adcSamples_me_gre;

    gz_blips_me_gre = prepare_gz_blips(params_me_gre);

    for i_interl = 1:params_me_gre.spi.interl
%         max_area = [gy_me_gre.area];
%         max_area = max(max_area);
%         duration = round(max_area/params.gen.lims.maxGrad/params.gen.lims.gradRasterTime)*params.gen.lims.gradRasterTime;
        duration = 5e-4;    % Fized to 0.5ms
        gxReph_me_gre(i_interl) = mr.makeTrapezoid('x',params_me_gre.gen.lims,'Area',gx_me_gre(1,i_interl).area.*-1,'Duration',duration);
        gyReph_me_gre(i_interl) = mr.makeTrapezoid('y',params_me_gre.gen.lims,'Area',gy_me_gre(1,i_interl).area.*-1,'Duration',duration);
    end

    % Flip angle, 3.5e-3 (approx rf + rephasing time)
    tr_tmp_me_gre = 3.5e-3+mr.calcDuration(gz_blips(1))+mr.calcDuration(mr.calcDuration(gx_me_gre)*params.gen.me_gre_echos) ... 
         +mr.calcDuration(gx_spoil)+mr.calcDuration(gxReph_me_gre); % aprox TR
    me_gre_tr_delay = params.gen.me_gre_tr-tr_tmp_me_gre;
    tr_tmp_me_gre = params.gen.me_gre_tr;
%     if and(params.gen.ro_type=='s',or(params.spi.type == 1,params.spi.type == 2)); tr_tmp_me_gre=tr_tmp_me_gre+mr.calcDuration(gx_pre_me_gre(1));end 
    params_me_gre = prepare_flip_angle(tr_tmp_me_gre,gx_me_gre, params_me_gre);
   
    % Create RF and Gz
    [rf0_me_gre,gz_me_gre,gzReph_me_gre,params_me_gre] =  create_rf_gz_pulseq(params_me_gre,tr_tmp_me_gre);
end

%% Prepare Delays and triggers
if params.gen.ro_type == 's'
    if params.spi.type == 0 && params.gen.te > 0
        te_delay = mr.makeDelay(round(params.gen.te-((mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph)+mr.calcDuration(gz_blips)),4));       % TE delay
        if params.gen.fid_nav == 1; te_delay.delay = te_delay.delay-(fid_nav_adc.duration*2)-mr.calcDuration(fid_nav_delay1)-mr.calcDuration(fid_nav_delay2); end
    elseif params.spi.type == 1 && params.gen.te > 0
        if params.gen.te > 0; te_delay = mr.makeDelay(round(params.gen.te-((mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph) ...
                +mr.calcDuration(gz_blips)+mr.calcDuration(gx(1))+mr.calcDuration(gx_pre(1))),4)); end
        if params.gen.fid_nav == 1; te_delay.delay = te_delay.delay-(fid_nav_adc.duration*2)-mr.calcDuration(fid_nav_delay1)-mr.calcDuration(fid_nav_delay2); end
    end
elseif params.gen.ro_type == 'c'
    if params.gen.te > 0; te_delay = mr.makeDelay(round(params.gen.te-(((mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph) ...
            +mr.calcDuration(gy_pre))+(mr.calcDuration(gx(1))*(ceil(params.epi.n_lines*params.epi.pf/params.epi.seg/2)+3))),4)); end         % TE delay +3 navig
end
if params.gen.tr_delay > 0; tr_delay = mr.makeDelay(params.gen.tr_delay); end                       % TR delay
if params.vaso.f_v_delay > 0; f_v_delay = mr.makeDelay(params.vaso.f_v_delay);          end         % FOCI-VASO delay
if params.vaso.v_b_delay > 0; v_b_delay = mr.makeDelay(params.vaso.v_b_delay);          end         % VASO-BOLD delay
if params.vaso.b_f_delay > 0; b_f_delay = mr.makeDelay(params.vaso.b_f_delay);          end         % BOLD-FOCI delay
if params.gen.me_gre > 0; me_gre_delay = mr.makeDelay(1); end                                       % Delay after ME GRE scan
if params.gen.me_gre > 0; me_gre_tr_delay = mr.makeDelay(me_gre_tr_delay); end                      % TR delay ME GRE scan
% TR delay for Cartesian, to avoid PNS
if params.gen.ro_type == 'c'
    if params.epi.tr_delay > 0
        tr_delay_epi = mr.makeDelay(params.epi.tr_delay);  
    end
end

% Skope delays
if params.gen.skope ~= 0
    % ToDo: Need to confirm the times of this delays
    sk_pre_delay    = mr.makeDelay(0.5);
%     sk_int_delay    = mr.makeDelay(400e-3);  % original
    sk_int_delay = mr.makeDelay(mr.calcDuration(gz_blips(1))-10e-6); % To compensate for missing gz blip in center partition
    sk_post_delay   = mr.makeDelay(3);
    sk_min_tr_delay = gzReph(1).riseTime+gzReph(1).flatTime+gzReph(1).fallTime + ...
                      gz_blips(1).riseTime + gz_blips(1).flatTime + gz_blips(1).fallTime ...
                      + mr.calcDuration(gx(1)); %+ mr.calcDuration(gx_ramp(1));
    sk_min_tr_delay = mr.makeDelay(120e-3-sk_min_tr_delay);                         % Here I take Skope minTR = 110ms
    % Delay to account for the missing RF pulses...
    sk_no_rf_delay = mr.calcDuration(gz_fs)+mr.calcDuration(rf0(1))+mr.calcDuration(gzReph(1));
    sk_no_rf_delay = mr.makeDelay(sk_no_rf_delay);
    if params.gen.field_strength == 7i
        skope_trig = mr.makeDigitalOutputPulse('ext1','duration',10e-6,'system',params.gen.lims);  % Skope trigger
    else
        skope_trig = mr.makeDigitalOutputPulse('osc0','duration',10e-6,'system',params.gen.lims);  % Skope trigger
    end   
    
end

% Multi Echo scan delays
if params.gen.seq == 3
    if params.gen.ro_type == 'c'
        tmp = (mr.calcDuration(rf0(1))/2)+mr.calcDuration(gzReph(1))+mr.calcDuration(gz_blips(1))+mr.calcDuration(gx_pre); % original
        tmp = tmp + (round(round(params.gen.n(2)/params.epi.ry)/2)*(mr.calcDuration(gx(1)))+mr.calcDuration(gy_blips(1)));
    elseif params.gen.ro_type == 's'
        if params.spi.type == 0
            tmp = round(((mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph)+mr.calcDuration(gz_blips)),4);       % TE delay
            if params.gen.fid_nav == 1; tmp = tmp+(fid_nav_adc.duration*2)+mr.calcDuration(fid_nav_delay1)+mr.calcDuration(fid_nav_delay2); end
        elseif params.spi.type == 1
%                 if params.gen.te > 0; te_delay = mr.makeDelay(round(params.gen.te-((mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph) ...
%                         +mr.calcDuration(gz_blips)+mr.calcDuration(gx(1))+mr.calcDuration(gx_pre(1))),4)); end
%                 if params.gen.fid_nav == 1; te_delay.delay = te_delay.delay-(fid_nav_adc.duration*2)-mr.calcDuration(fid_nav_delay1)-mr.calcDuration(fid_nav_delay2); end
        end
    end
    for i=1:params.gen.echos
        fm_te_delay0(i) = mr.makeDelay(round(params.gen.multi_te(i)-tmp,4));
    end
    for i=1:params.gen.echos-1
        fm_te_delay1(i) = mr.makeDelay(round(fm_te_delay0(1).delay-fm_te_delay0(i+1).delay,4));
    end
    if params.gen.seq == 3; multi_te_delay = mr.makeDelay(params.gen.multi_te(1)-params.gen.multi_te(2)); end   % Multi TE delay
end

% Delay for repetitions without MT pulse, only for ABC (seq=2)
if params.gen.seq == 2
    no_mt_delay = mr.makeDelay(mr.calcDuration(MT)+mr.calcDuration(gx_spoil));
end

% No Gz blip delay
no_blip_delay = mr.makeDelay(mr.calcDuration(gz_blips(1)));

% No FatSat delay...
if params.gen.fat_sat
    no_fs_delay = mr.makeDelay(mr.calcDuration(rf_fs));
end

% External trigger for fMRI
ext_trig = mr.makeDigitalOutputPulse('ext1','duration',10e-6,'system',params.gen.lims);        % External trigger
dummy_delay = mr.makeDelay(3000e-3);                                                  % Delay to use as dummy anywhere

% Save partitions and segments in a new variable
last_part = params.gen.n(3);
last_seg = params.gen.seg;
if params.gen.seq == 3
    last_contrast = params.gen.echos;
    last_ech = 1;
else
    last_contrast = 1;
    last_ech = params.gen.echos;
end

%% Add blocks to Skope seq
if params.gen.skope == 1 || params.gen.skope == 3
    seq_sk=mr.Sequence();          % Create a new sequence object
%     seq_sk.addBlock(sk_pre_delay); % Skope initial delay
    % Skope sync scans
    for i_sk_pre = 1:params.gen.skope_sync
        seq_sk.addBlock(mr.makeLabel('SET','ONCE',1));
%         seq_sk.addBlock(sk_pre_delay);
        seq_sk.addBlock(skope_trig);
        seq_sk.addBlock(sk_int_delay);     % Gradient free interval
        seq_sk.addBlock(no_blip_delay);
        if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre(1)));
            seq_sk.addBlock(tmp_delay);                       
        end
        tmp_delay = mr.makeDelay(mr.calcDuration(gx(1)));
        seq_sk.addBlock(tmp_delay,adc);
        tmp_delay = mr.makeDelay(mr.calcDuration(gz_spoil));
        seq_sk.addBlock(tmp_delay);
        % rewinder if multiple echos
        if params.gen.echos > 1
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre(1)));
            seq_sk.addBlock(tmp_delay)
        end
        % Delay after sync scans
        if i_sk_pre == params.gen.skope_sync; seq_sk.addBlock(sk_post_delay); end
        seq_sk.addBlock(mr.makeLabel('SET','ONCE',0));
    end
    % Partitions loop
    for i=1:last_part
        % Spiral/Cartesian Interleaves/segments loop
        % for m=1:params.gen.n(3)
            for j = 1:last_seg
                adc.phaseOffset = adc_phase_offset(i);
%                 seq_sk.addBlock(gzReph(i));
                if params.gen.te > 0; seq_sk.addBlock(te_delay); end % TE delay
                % EPI navigators
                if params.gen.ro_type == 'c'
                    for i_nav = 1:3
                        gx.amplitude = -gx.amplitude;
                        seq_sk.addBlock(gx,adc); 
                    end
                end
                seq_sk.addBlock(sk_no_rf_delay);      % Delay to account for the missing RF pulses... 
%                 Trigger Before readout...
%                 We only add the trigger in the center partition, remove if if I want it in all partitions
                if i == floor((params.gen.n(3)/2)+1)
                    sk0 = seq_sk.duration();           % To get Skope triger-ADC delay
                    % Trigger before spiral readout..
                    seq_sk.addBlock(skope_trig);
                    seq_sk.addBlock(sk_int_delay);     % Gradient free interval
%                     seq_sk.addBlock(dummy_delay);      % Temp: dummy delay to record oscilations after readout...
                end
                if params.gen.ro_type == 's'
                    seq_sk.addBlock(gz_blips(i,j));
                    % if spiral-in or in-out, we need in-plane pre-phasing
                    if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
                        seq_sk.addBlock(gx_pre(i,j),gy_pre(i,j));                       
                    end
                    if i == floor((params.gen.n(3)/2)+1); sk1 = seq_sk.duration();    end                % To get Skope triger-ADC delay
                    seq_sk.addBlock(gx(i,j),gy(i,j),adc);
                    seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
%                     %%%% Temp: Trigger after spiral readout
%                     if i == floor((params.gen.n(3)/2)+1)
%                         seq_sk.addBlock(skope_trig);
%                         seq_sk.addBlock(sk_int_delay);     % Gradient free interval
%                     end
% %                     seq_sk.addBlock(dummy_delay);      % Temp: dummy delay to record oscilations after readout...
%                     %%%%
                    % rewinder if multiple echos
                    if params.gen.echos > 1; seq_sk.addBlock(gx_pre(i,j),gy_pre(i,j)); end
                    if params.gen.tr_delay > 0; seq_sk.addBlock(tr_delay); end % TR delay
                elseif params.gen.ro_type == 'c'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq_sk.addBlock(gx_pre,gy_pre(j))
                    else
                        seq_sk.addBlock(gx_pre,gy_pre(j),gz_blips(i,j))
                    end
                    for k = 1:round(params.gen.n(2)/params.epi.ry)-((round(params.gen.n(2)/params.epi.ry)-round(params.gen.n(2)/params.epi.ry*params.epi.pf))/2)
                        gx.amplitude = -gx.amplitude;
                        if k == 1
                            seq_sk.addBlock(gx,gy_blip_up,adc);
                        elseif k == round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
                            seq_sk.addBlock(gx,gy_blip_down,adc);
                        else
                            seq_sk.addBlock(gx,gy_blips,adc);
                        end
                    end
                    if gx.amplitude > 0
                        gx.amplitude = -gx.amplitude;
                    end
                    % seq_sk.addBlock(gx_spoil,gy_spoil,gz_spoil);
                    if params.gen.tr_delay > 0; seq_sk.addBlock(tr_delay); end % TR delay
                    if params.epi.tr_delay > 0; seq_sk.addBlock(tr_delay_epi); end
                end
%                 seq_sk.addBlock(sk_min_tr_delay); % Do I need this if I only trigger center partition?
            end
        % end
    end
    params.gen.skope_trig_adc_delay = sk1-sk0;
end

%% Actual scan
seq = mr.Sequence();
% ME-GRE calibration
if params.gen.me_gre > 0
    if params.gen.me_gre == 1; seq.addBlock(mr.makeLabel('SET','ONCE',1)); end
%     % FOCI
%     if params.gen.seq == 1
%         if params.vaso.foci; seq.addBlock(rf_foci); end              % FOCI
%         if params.vaso.f_v_delay > 0; seq.addBlock(f_v_delay); end   % FOCI-VASO delay
%     end
    for i=1:params_me_gre.gen.n(3)
        for j=1:params_me_gre.spi.interl
            % Fat Sat
            if params_me_gre.gen.fat_sat       % fat-sat
                seq.addBlock(gx_fs,gy_fs,gz_fs);    % spoilers
                seq.addBlock(rf_fs);
%                 gx_fs.amplitude = gx_fs.amplitude*-1; gy_fs.amplitude = gy_fs.amplitude*-1; gz_fs.amplitude = gz_fs.amplitude*-1;
%                 seq.addBlock(gx_fs,gy_fs,gz_fs);    % spoilers
            end
            % Spoil
            seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
            rf = rf0_me_gre(i);
            rf.phaseOffset = rf_phase_offset_me_gre(i,j);
            adc_me_gre.phaseOffset = adc_phase_offset_me_gre(i,j);
            adc_post.phaseOffset = adc_phase_offset_me_gre(i,j);
            if i==1 && j==1; te0_me_gre = seq.duration(); end         % save seq dur to calc TE
            if i==1 && j==1; tr0_me_gre = seq.duration(); end         % save seq dur to calc TR
            seq.addBlock(rf,gz_me_gre(i));
            seq.addBlock(gzReph_me_gre(i));
            if params_me_gre.spi.type == 1 || params_me_gre.spi.type == 3 || params_me_gre.spi.type == 4
                seq.addBlock(gx_pre(i,j),gy_pre(i,j));                       
            end
            seq.addBlock(gz_blips_me_gre(i,j));
            for k=1:params_me_gre.gen.echos
                % Save TEs
                if i==1 && j==1; params_me_gre.gen.te(k) = seq.duration()-te0_me_gre-(mr.calcDuration(rf0_me_gre)/2); end
                seq.addBlock(gx_me_gre(i,j),gy_me_gre(i,j),adc_me_gre);
                seq.addBlock(gxReph_me_gre(j),gyReph_me_gre(j));
%                 if k==params_me_gre.gen.echos; seq.addBlock(gx_spoil,gy_spoil,gz_spoil);end
            end
            seq.addBlock(me_gre_tr_delay); % TR delay
            if i==1 && j==1; tr1_me_gre = seq.duration(); end                           % save seq dur to calc TR
        end
    end
    seq.addBlock(me_gre_delay); % Delay after the GRE ME scan..
    if params.gen.me_gre == 1; seq.addBlock(mr.makeLabel('SET','ONCE',0)); end
end
% If ME GRE separate scan, split the sequence
if params.gen.me_gre == 2
    seq_me_gre = seq;
    seq = mr.Sequence();
end

% Skope
if params.gen.skope == 2 || params.gen.skope == 3
    % Skope sync scans
    for i_sk_pre = 1:params.gen.skope_sync
        seq.addBlock(mr.makeLabel('SET','ONCE',1));
        seq.addBlock(sk_pre_delay);
        seq.addBlock(skope_trig);
        seq.addBlock(sk_int_delay);     % Gradient free interval
        seq.addBlock(no_blip_delay);
        if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre));
            seq.addBlock(tmp_delay);                       
        end
        tmp_delay = mr.makeDelay(mr.calcDuration(gx(1)));
        seq.addBlock(tmp_delay,adc);
        tmp_delay = mr.makeDelay(mr.calcDuration(gz_spoil));
        seq.addBlock(tmp_delay);
        % rewinder if multiple echos
        if params.gen.echos > 1
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre));
            seq.addBlock(tmp_delay) 
        end
        % Delay after sync scans
        if i_sk_pre == params.gen.skope_sync; seq.addBlock(sk_post_delay); end
        seq.addBlock(mr.makeLabel('SET','ONCE',0));
    end
end

seq_t0 = seq.duration();
% MT specific 
if params.gen.seq == 2
    if params.mt.mt_prep > 0 && params.mt.mt == 1
%         seq.addBlock(mr.makeLabel('SET','ONCE',1));
        for i_mt_prep = 1:params.mt.mt_prep
            seq.addBlock(MT);                                      % MT pulse 
            seq.addBlock(gz_spoil);
        end
        seq.addBlock(mr.makeDelay(1));                            % Delay after dummy MT-prep 
%         seq.addBlock(mr.makeLabel('SET','ONCE',0));
    end
    if params.mt.rf_spoil == 0
        rf_phase_offset = rf_phase_offset.*0;
        adc_phase_offset = adc_phase_offset.*0;
    end
end
% seq.addBlock(mr.makeLabel('SET','REP',0));
% FOCI
seq.addBlock(ext_trig);                                      % External trigger
if params.gen.seq == 1
    if params.vaso.foci; seq.addBlock(rf_foci); end              % FOCI
    if params.vaso.f_v_delay > 0; seq.addBlock(f_v_delay); end   % FOCI-VASO delay
end

% Main loop
for i_ro_blocks = 1:ro_blocks
    for i=1:last_part
        % Adding MT pulse every specified time
        if params.gen.seq == 2
            if params.mt.mt == 1 
                    seq.addBlock(MT);                                      % MT pulse 
                    seq.addBlock(gz_spoil);                                 % Original
            else
                    seq.addBlock(no_mt_delay);
            end
        end
        for l=1:last_contrast
            for j=1:last_seg
                % Fat Sat
                if (params.gen.fat_sat && params.gen.fs_interl==1) || (params.gen.fat_sat && params.gen.fs_interl==0 && j==1)
                    seq.addBlock(gx_fs,gy_fs,gz_fs);    % spoilers
                    seq.addBlock(rf_fs);
                elseif params.gen.fat_sat && params.gen.fs_interl==0 && j>1
                    seq.addBlock(gx_fs,gy_fs,gz_fs);    % spoilers
                    seq.addBlock(no_fs_delay); 
                end
                % Spoil
                seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                if i==1 && j==1 && i_ro_blocks==1; tr0 = seq.duration(); end         % save seq dur to calc TR
                rf = rf0(i);
                rf.phaseOffset = rf_phase_offset(i,j);
                adc.phaseOffset = adc_phase_offset(i,j);
                adc_post.phaseOffset = adc_phase_offset(i,j);
                seq.addBlock(rf,gz(i));
                if i==1 && j==1; te0 = seq.duration(); end                           % save seq dur to calc TE
                seq.addBlock(gzReph(i));
                % EPI navigators
                if params.gen.ro_type == 'c'
                    if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
                    for i_nav = 1:3
                        gx.amplitude = -gx.amplitude;
                        seq.addBlock(gx,adc); 
                    end
                end
                % FID navigators
                if params.gen.ro_type == 's' && params.gen.fid_nav
                    if i==1 && j==1; te0_fid = seq.duration(); end
                    seq.addBlock(fid_nav_delay1);
                    if i==1 && j==1; te1_fid = seq.duration(); end
                    seq.addBlock(fid_nav_adc);
                    seq.addBlock(fid_nav_delay2);
                    if i==1 && j==1; te2_fid = seq.duration(); end
                    seq.addBlock(fid_nav_adc);
                    params.gen.fid_nav_te = [te1_fid-te0_fid,te2_fid-te0_fid]+(mr.calcDuration(rf0)/2)+(mr.calcDuration(fid_nav_adc)/2);
                end
                % Echos loop
                for k=1:last_ech
                    if params.gen.ro_type == 's'
                        if k==1 || k==params.gen.echos-1
                            if i == 1 && params.gen.kz_enc == 1
                                if params.gen.skope == 2 || params.gen.skope == 3
                                    seq.addBlock(skope_trig);
                                    seq.addBlock(sk_int_delay);     % Gradient free interval
                                else
                                    seq.addBlock(no_blip_delay);
                                    if i_ro_blocks==1; params.vaso.ti1 = seq.duration();end
                                    if i_ro_blocks==2; params.vaso.ti2 = seq.duration();end
                                end
                            else
                                seq.addBlock(gz_blips(i,j));
                            end
                        end
                        if and(params.gen.seq == 3, params.gen.multi_te(1) > 0); seq.addBlock(fm_te_delay0(l)); end
                        if and(params.gen.seq ~= 3, params.gen.te > 0); seq.addBlock(te_delay); end       % TE delay
                        % if spiral-in or in-out, we need in-plane pre-phasing
                        if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
                            seq.addBlock(gx_pre(i,j),gy_pre(i,j));                       
                        end
                        if i==1 && params.spi.type==0 && j==1; te1 = seq.duration(); end                           % save seq dur to calc TE
                        seq.addBlock(gx(i,j),gy(i,j),adc);
    %                     seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                        if i==1 && params.spi.type==1 && j==1; te1 = seq.duration(); end                           % save seq dur to calc TE
                        if i==1 && (params.spi.type==3 || params.spi.type==4) && j==1; te1 = seq.duration()-(mr.calcDuration(gx(1)))/2; end                           % save seq dur to calc TE
                        if params.gen.dork; seq.addBlock(adc_post); end
%                         % Spoil
%                         if params.spi.interl > 1  % Not spoiling in all direction with Multi-shot??
% %                             if k==params.gen.echos; seq.addBlock(gz_spoil); end  % Original
%                             seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
%                         else
% %                             if k==params.gen.echos; seq.addBlock(gz_spoil,gx_spoil);end  % Original
%                             seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
%                         end
                        % rewinder if multiple echos
                        if params.gen.echos > 1 && params.gen.seq == 4; seq.addBlock(gx_pre(i,j),gy_pre(i,j)); end
                        if params.gen.tr_delay > 0; seq.addBlock(tr_delay); end % TR delay
                        if (params.gen.seq == 3 && params.gen.multi_te(1) > 0 && l > 1)
                            seq.addBlock(fm_te_delay1(l-1));
                        end
                        %%% I want to put the second delay here....
                    elseif params.gen.ro_type == 'c'
                        if i == floor((params.gen.n(3)/2)+1)
                            seq.addBlock(gx_pre,gy_pre(j))
                        else
                            seq.addBlock(gx_pre,gy_pre(j),gz_blips(i,j))
                        end
                        for k = 1:params.epi.n_lines/params.epi.seg
                            gx.amplitude = -gx.amplitude;
                            if k == 1
                                seq.addBlock(gx,gy_blip_up,adc);
                            elseif k == params.epi.n_lines/params.epi.seg
                                seq.addBlock(gx,gy_blip_down,adc);
                            else
                                seq.addBlock(gx,gy_blips,adc);
                            end
                            if i==1 && j==1 && k==ceil((params.epi.n_lines/params.epi.seg*params.epi.pf/2)); te1=seq.duration(); end                           % save seq dur to calc TE
                        end
                        if gx.amplitude > 0
                            gx.amplitude = -gx.amplitude;
                        end
                        seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                        if params.gen.tr_delay > 0; seq.addBlock(tr_delay); end % TR delay
                        if params.epi.tr_delay > 0; seq.addBlock(tr_delay_epi); end
                    end
                end
                if i==1 && j==1 && i_ro_blocks==1; tr1 = seq.duration(); end                           % save seq dur to calc TR
            end
        end
%         seq.addBlock(mr.makeLabel('INC','PAR',1));
        if params.gen.seq == 3 && i_ro_blocks > 1
%             seq.addBlock(multi_te_delay); 
        end                                % Multi-TE delay to match TRs
        dur1 = seq.duration();
    end
    if params.gen.seq == 1 && i_ro_blocks == 1
        if params.vaso.v_b_delay > 0; seq.addBlock(v_b_delay); end      % VASO-BOLD delay
    elseif params.gen.seq == 1 &&  i_ro_blocks == 2
        if params.vaso.b_f_delay > 0; seq.addBlock(b_f_delay); end   % BOLD-FOCI delay
    end
% seq.addBlock(dummy_delay);
%     seq.addBlock(mr.makeLabel('INC','REP',1));
%     seq.addBlock(mr.makeLabel('SET','PAR',0));
end
seq_t1 = seq.duration();

%% Check timing of the sequence
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Getting k-space trajectories
traj_st = 1;
if params.gen.me_gre == 1
    ks_traj_me_gre = create_ks_trajectory(seq,adc_me_gre,params_me_gre,traj_st,0);
    traj_st = params_me_gre.gen.ro_samples*params_me_gre.spi.interl*params_me_gre.gen.echos*params_me_gre.gen.n(3)+1;
elseif params.gen.me_gre == 2
    ks_traj_me_gre = create_ks_trajectory(seq_me_gre,adc_me_gre,params_me_gre,traj_st,0);
    traj_st = 1;
end
if params.gen.ro_type == 'c'; traj_st = (adc.numSamples*3)+1; end
if params.gen.ro_type == 's' && params.gen.fid_nav == 0; traj_st = 1; end
if params.gen.ro_type == 's' && params.gen.fid_nav == 1; traj_st = 201; end % The FIDs are 100samples each, we have two
ks_traj = create_ks_trajectory(seq,adc,params,traj_st,1);

%% Adding some extra parameters to params directory
if params.gen.me_gre > 0
    params_me_gre = prepare_add_parameters(ks_traj_me_gre,gx_me_gre,rf,adc_me_gre,te0_me_gre,te1,tr0_me_gre,tr1_me_gre,seq_t0,seq_t1,params_me_gre);
end
params = prepare_add_parameters(ks_traj,gx,rf,adc,te0,te1,tr0,tr1,seq_t0,seq_t1,params);

%% Check accoustic resonance frequencies, Des script
check_accoustic_fq_pns(seq,params,seq_t0)

%% Set definitions
seq.setDefinition('MaxAdcSegmentLength',params.gen.adc_split);
seq.setDefinition('FOV', params.gen.fov);
seq.setDefinition('Name', seq_name);
seq.setDefinition('baseResolution', 2000);
if params.gen.me_gre == 2
    seq_me_gre.setDefinition('MaxAdcSegmentLength',params_me_gre.gen.adc_split);
    seq_me_gre.setDefinition('FOV', params_me_gre.gen.fov);
    seq_me_gre.setDefinition('Name', sprintf('%s_me_gre',seq_name));
    seq_me_gre.setDefinition('baseResolution', 64);
end

%% Saving files
% Check if folder exist, if no it creates it
if exist('./data','dir') == 0; mkdir ./data/; end
if exist(sprintf('./data/%s',folder_name),'dir') == 0
    tmp = sprintf('./data/%s',folder_name);
    system(sprintf('mkdir %s',tmp));
    system(sprintf('mkdir %s/acq',tmp));
    system(sprintf('mkdir %s/analysis',tmp));
    system(sprintf('mkdir %s/ismrmd',tmp));
    system(sprintf('mkdir %s/ismrmd/2d',tmp));
    system(sprintf('mkdir %s/ismrmd/3d',tmp));
    system(sprintf('mkdir %s/raw',tmp));
    system(sprintf('mkdir %s/raw/twix',tmp));
    system(sprintf('mkdir %s/raw/dicom',tmp));
    system(sprintf('mkdir %s/raw/nifti',tmp));
    system(sprintf('mkdir %s/recon',tmp));
    system(sprintf('mkdir %s/recon/2d',tmp));
    system(sprintf('mkdir %s/recon/3d',tmp));
    system(sprintf('mkdir %s/tmp',tmp));
    system(sprintf('mkdir %s/sim',tmp));
end

namestr = strcat('./data/',folder_name,'/acq/',seq_name);
seq.write(strcat(namestr,'.seq'));
if params.gen.skope == 1 || params.gen.skope == 3
    seq_sk.write(strcat(namestr,'_sk.seq'));
end
if params.gen.me_gre == 2
    seq_me_gre.write(strcat(namestr,'_me_gre.seq'));
end
if params.spi.type == 3 && params.gen.ro_type == 's'
    % ToDo : Fix this part so it works for more echos..
    ks_traj_full = ks_traj;
    ks_traj = ks_traj_full.e1;
    save(strcat(namestr,'_e1_ks_traj_nom.mat'),'ks_traj')
    ks_traj = ks_traj_full.e2;
    save(strcat(namestr,'_e2_ks_traj_nom.mat'),'ks_traj')
else
    save(strcat(namestr,'_ks_traj_nom.mat'),'ks_traj')
end
if params.gen.me_gre > 0
    save(strcat(namestr,'_ks_traj_me_gre_nom.mat'),'ks_traj_me_gre')
    save(strcat(namestr,'_params_me_gre.mat'),'params_me_gre')
end
save(strcat(namestr,'_params.mat'),'params')
% seq.plot()
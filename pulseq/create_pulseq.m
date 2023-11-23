clear all; clc
% Change directory to path where sosp_vaso git has been cloned
cd /home/amonreal/Documents/PhD/PhD_2023/sosp_vaso/
addpath(genpath("./pulseq/functions"))
% Add path to pulseq installation folder
addpath(genpath("/home/amonreal/Documents/PhD/tools/pulseq_1.4.0/"))
addpath(genpath("/home/amonreal/Documents/PhD/tools/tOptGrad_V0.2/minTimeGradient/mex_interface/"))
addpath(genpath("/home/amonreal/Documents/PhD/tools/pns_prediction/"))
addpath(genpath("/home/amonreal/Documents/PhD/tools/check_grad_idea_Des/"))
warning('OFF', 'mr:restoreShape')

%% ToDo
% - Check exactly how the saturation pulse should be, what spoliers I need?
% - How am I gonna take into account the gradient delays (trajectory,adc
% and time vector)
% - Check calculation of TE,TR and all params for all sequences
% - Do I need to include adc delay? as they do in some pulseq examples

%% Define parameters
folder_name = '11212023_abc';        % Day I am scanning
seq_name = 'abc_02';                % use sv/abc/sb_n (n for the diff scans at each day).. sb_04_OUT_08mm_2seg
params.gen.seq = 4;                 % 1-VASO 2-ABC 3-Multi-Echo 4- BOLD
params.gen.field_strength = 7;      % Field Strength (7=7T,7i=7T-impuse_grad,9=9.4T,11=11.7T)
params.gen.traj_scan = 0;           % Generate sequence to measure trajectory.

% General parameters
params.gen.fov = [192 192 24].*1e-3; % (192,192,24), reduced FOV=(140,140,24)
params.gen.res = [0.8 0.8 1].*1e-3; % Target resoultion 
params.gen.fa = 0;                  % Set to 0, to use Ernst Angle
params.gen.ernst_t1 = 1800e-3;      % T1 to calculate Ernst Angle (s) 7T=(1220e-3)(1800e-3) 9T=(1425e-3)(2100e-3)
params.gen.te = 10e-3;               % Set to 0 to shortest TE possible
params.gen.tr_delay = 0e-3;         % Delay between acquisitions in sec
params.gen.ro_type = 's';           % 's'-Spiral, 'c'-Cartesiaen
params.gen.kz = 1;                  % Acceleration in Kz
params.gen.pf = 1;                  % Partial fourier in Kz
params.gen.fat_sat = 1;             % Fat saturation (1=yes,0=no)
params.gen.vfa = 0;                 % Variable FA, demo
params.gen.skope = 0;               % Add skope sync scan and triggers, 0=N0, 1=sepscan, 2=same scan/concurrent(center partition)
params.gen.dork = 0;                % extra adc's for DORK correction
params.gen.kz_enc = 1;              % k-space partition encoding 0=linear,1=center-out For Cartesian now only linear encoding
params.gen.ph_oversampling = 0;     % Partition phase oversampling in %, to avoid partition phase-wrap (10)
params.gen.echos = 1;               % Echos per RF pulse
         
% Spiral parameters
params.spi.type = 0;                % spiral type 0=spiral-Out , 1=spiral-In, 3=In-Out
params.spi.rotate = 'none';         % Spiral rotation ('none','golden','180','120','linear'), linear not implemented
params.spi.increment = 'linear';    % Spiral increment mode (for now only linear)
params.spi.max_grad  = 40;          % Peak gradient amplitude for spiral (mT/m)  (7T=35) (9T=50) (7i=75) (7T/6int=40)
params.spi.max_sr = 155;            % Max gradient slew rate for spiral (mT/m/ms) (7T=155) (9T=250) (7i=750) (7T/6int=155)
params.spi.interl = 6;              % Spiral interleaves
params.spi.vd = 1.6;                % Variability density
params.spi.rxy = 3;                 % In-plane undersampling
params.spi.bw = 500e3;              % Spiral BW in Hz (Max value 1,000e3) (500e3)

% MT pulse parameters
params.mt.mt = 1;                   % Add MT pulse, 0 for reference scan without MT
params.mt.alpha = 225;  %225
params.mt.delta = -650;  %650      % for Pulseq approach should be 650 to match Viktors phase, normally -650
params.mt.trf = 0.004;
params.mt.mt_rep = 100e-3;          % How often to play the mt pulse (130e-3)
params.mt.bold = 0;                 % Get BOLD reference acq after mt one

% EPI parameters
params.epi.ry = 4;
params.epi.pf = 1;                % In-plane PF, (1,7/8,6/8)
params.epi.te = [33.6 36.6 38.6 40]*1e-3+0.07; % Echo times for field map
params.epi.seg = 1;                 % EPI Segments
params.epi.tr_delay = 0;            % Delay after each TR, needed to reduce SAR
params.epi.bw_px = 1046;             % BW in Hz rxy=3,pf=6/8->960 (1046)

% VASO parameters
params.vaso.foci = 1;               % FOCI inversion?
params.vaso.bold_ref = 1;           % BOLD reference volume
params.vaso.tr = 4500e-3*2;         % volume TR (4500e-3)
params.vaso.ti1 = 1800e-3;          % VASO TI1, 
params.vaso.ti2 = params.vaso.ti1+(params.vaso.tr/2);
params.vaso.foci_ampl = 270;        % FOCI amplitude (140/270)
params.vaso.foci_dur = 10410e-6;    % FOCI duration
params.vaso.foci_bw = 150;          % FOCI Bandwidth (150)
params.vaso.f_v_delay = 600e-3;     % FOCI-VASO delay (600e-3)
params.vaso.v_b_delay = 10e-3;      % VASO-BOLD delay (10e-3)
params.vaso.b_f_delay = 5e-3;       % BOLD-FOCI delay (5e-3)

% Set system limits
if params.gen.field_strength == 7
    % 7T
    lims = mr.opts('MaxGrad',65,'GradUnit','mT/m',...
        'MaxSlew',200,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 100e-6,'adcDeadTime', 10e-6, 'B0',6.98);  % To read it in VM I need rfDeadtime = 180e-6, 100e-6 for scanner
elseif params.gen.field_strength == 9
    % 9.4T
    lims = mr.opts('MaxGrad',80,'GradUnit','mT/m',...
        'MaxSlew',333,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 150e-6,'adcDeadTime', 10e-6, 'B0',9.38);  % To read it in VM I need rfDeadtime = 180e-6
elseif params.gen.field_strength == 7i
    % 7T impuse gradient
    lims = mr.opts('MaxGrad',198,'GradUnit','mT/m',...
        'MaxSlew',900,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 100e-6,'adcDeadTime', 10e-6, 'B0',6.98);  % To read it in VM I need rfDeadtime = 180e-6
end

% Set gradient files
if params.gen.field_strength == 7
    grad_file = '/home/amonreal/Documents/PhD/IDEA/gradient_files/MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc';
elseif params.gen.field_strength == 7i
    grad_file = '/home/amonreal/Documents/PhD/IDEA/gradient_files/MP_GPA_K2298_2250V_1250A_AC207_Base.asc';
elseif params.gen.field_strength == 9
    grad_file = '';
end

% Some calculations, restrictions in seq...
params.gen.del_k = (1./params.gen.fov);
if params.gen.ro_type == 's'; params.gen.del_k(1:2) = params.gen.del_k(1:2).*params.spi.rxy; end
if params.gen.ro_type == 'c'; params.spi.interl = 1; end
if params.gen.ro_type == 'c'; params.gen.del_k(2) = params.gen.del_k(2).*params.epi.ry; end
if params.mt.bold; params.mt.mt = 1; end        % If MT BOLD corrected, use MT 
params.gen.del_k(3) = params.gen.del_k(3).*params.gen.kz;
params.gen.n = round(params.gen.fov./params.gen.res);
if and(params.gen.seq == 1,params.vaso.bold_ref == 1) || and(params.gen.seq == 2,params.mt.bold)
    ro_blocks = 2; 
elseif params.gen.seq == 3
    ro_blocks = length(params.epi.te);
else
    ro_blocks = 1; 
end
% Partition or phase oversampling
params.gen.n_ov = params.gen.n;
if params.gen.ph_oversampling > 0
    params.gen.del_k(3) = params.gen.del_k(3)/(1+(params.gen.ph_oversampling/100)); 
    params.gen.n(3) = round(round(params.gen.n(3)*(1+(params.gen.ph_oversampling/100)))/2)*2;
end

% Trying to make n multiple of 4,update res
tmp = mod(params.gen.n,4);
params.gen.n(1:2) = params.gen.n(1:2)+tmp(1:2);
% ToDo: Check... Not sure if I do need this...
params.gen.n(3) = params.gen.n(3)/params.gen.kz/params.gen.pf;
% Making sure fovz/rz are integers and even
if round((params.gen.fov(3)/params.gen.kz)*1000,9)/round((params.gen.fov(3)/params.gen.kz)*1000) ~= 1
    params.gen.fov(3) = round(params.gen.fov(3)*1000/2/params.gen.kz)*2*params.gen.kz*1e-3;
end
% if params.gen.seq == 3; params.gen.ro_type = 'c'; end   % if fieldmap, cartesian
if params.gen.seq == 2; params.gen.ro_type = 's'; end   % if ABC, Spiral
if params.mt.mt == 0; params.mt.bold = 0; end           % if no MT, no ref BOLD  
if params.gen.ro_type == 'c'
    params.gen.seg = params.epi.seg;
elseif params.gen.ro_type == 's'
    params.gen.seg = params.spi.interl;
end


%% Create blocks (I will get this into functions)
% TR-FOCI
if params.gen.seq == 1
    [B1_foci,phase,rf_complex,Gz_foci,fa] = tr_foci(params);
    rf_foci = mr.makeArbitraryRf(rf_complex,fa*pi/180, 'system', lims);%, 'Delay',2e-4);
end

% MT pulse
if params.gen.seq == 2
%     % Using Viktor's code
%     vpulse = VPF_gaussian_pulse_4_Maastricht(params.mt.alpha,params.mt.delta,params.mt.trf);
%     MT = mr.makeArbitraryRf(vpulse.b1,params.mt.alpha*pi/180, 'system', lims);%, 'Delay',2e-4);
    % Using Pulseq
    MT = mr.makeGaussPulse(params.mt.alpha*pi/180,lims,'Duration',params.mt.trf,'FreqOffset',params.mt.delta);
end

% Fat sat:
if params.gen.fat_sat
    % AMM: ToDo: how to properly desing the fat-sat pulse
    sat_ppm = -3.45;
    if params.gen.field_strength == 7  || params.gen.field_strength == 7i 
        fs_angle = 80; % 110/80
    elseif params.gen.field_strength == 9
        fs_angle = 40; % 110/80
    end
    sat_freq = sat_ppm*1e-6*lims.B0*lims.gamma;
    if params.gen.field_strength == 7 || params.gen.field_strength == 7i 
        rf_fs = mr.makeGaussPulse(fs_angle*pi/180,'system',lims,'Duration',5e-3,...
            'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
    elseif params.gen.field_strength == 9
        rf_fs = mr.makeGaussPulse(fs_angle*pi/180,'system',lims,'Duration',8e-3,...
            'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
    end
    
    % AMM: ToDo: Check correct values for spoiling gradient
    gx_fs = mr.makeTrapezoid('x',lims,'MaxGrad',18e-3*lims.gamma,'Area',-1/1e-3,'Duration',1.5e-3);
    gy_fs = mr.makeTrapezoid('y',lims,'MaxGrad',18e-3*lims.gamma,'Area',1/1e-3,'Duration',1.5e-3);
    gz_fs = mr.makeTrapezoid('z',lims,'MaxGrad',18e-3*lims.gamma,'Area',1/1e-3,'Duration',1.5e-3);
    gx_fs.delay = mr.calcDuration(rf_fs);
    gy_fs.delay = mr.calcDuration(rf_fs);
    gz_fs.delay = mr.calcDuration(rf_fs);
end

% Kz Blips:
tmp = params.gen.n(3);
if mod(tmp,2) == 1; tmp = tmp + 1; end
for i=1:tmp
    if mod(params.gen.n(3),2) == 1
        area = -(params.gen.del_k(3)*(params.gen.n(3)/2))+(params.gen.del_k(3)*(i-1))+(params.gen.del_k(3)/2);
    else
        area = -(params.gen.del_k(3)*(params.gen.n(3)/2))+(params.gen.del_k(3)*(i-1));
    end
    dur = ceil(2*sqrt(area/lims.maxSlew)/10e-6)*10e-6;
    if area ~= 0
        % I fix the duration of the blip to an even number
%         gz_blips(i) = mr.makeTrapezoid('z',lims,'Area',area,'Duration',6e-4);
        gz_blips(i) = mr.makeTrapezoid('z',lims,'Area',area);
    end
end
% if mod(floor(params.gen.n(3)/params.gen.kz),2) == 1  % Original
if mod(floor(params.gen.n(3)),2) == 1
    gz_blips(round(params.gen.n(3)/2)) = [];    % Removing the empty blip...
else
    gz_blips(round(params.gen.n(3)/2)+1) = [];    % Removing the empty blip...
end

% Reshuffling blips if center-out
if params.gen.kz_enc == 1
    tmp = [];
    j = floor(params.gen.n(3)/params.gen.kz)/2;
    for i=1:length(gz_blips)
        tmp = [tmp gz_blips(j)];
        if mod(i,2) == 0
            j = j-i;
        else
            j = j+i;
        end
    end
    gz_blips = tmp;
end

% Spoilers
% AMM: Todo: Need to confirm the size (amplitude/area) of this spoiler
mag_spoil = 40e-3;     % mT
sr_spoil = 180;        % mT/m/s
gx_spoil=mr.makeTrapezoid('x','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);%,'delay',1e-3);
gy_spoil=mr.makeTrapezoid('y','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);%,'delay',1e-3);
% gz_spoil=mr.makeTrapezoid('z','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);%,'delay',1e-3); % original
gz_spoil=mr.makeTrapezoid('z',lims,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);%,'delay',1e-3);

%% Preparing readout elements
[rf_phase_offset,adc_phase_offset] = rf_adc_phase(params);
if params.gen.ro_type == 's'
    [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc1(params,lims);
        for j=1:params.gen.n(3)
            for i=1:params.spi.interl
        
                % Phase encoding delay
                if params.gen.kz_enc == 0
                    delay = mr.calcDuration(gz_blips(1)); 
                elseif params.gen.kz_enc == 1
                    delay = mr.calcDuration(gz_blips(end)); 
                end
%                 if j == 1 && i == 1 && params.spi.type == 0
%                     % Temp, removing sampled added to round gradient
%                     adcSamples = adcSamples-(delay/adcDwell);
%                 end

                % Readout gradients
                % Temp: for Pulseq version 1.4.1
%                     gx(j,i) = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:,i,j), lims, 'first', 0, 'last', 0);
%                     gy(j,i) = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:,i,j), lims, 'first', 0, 'last', 0);
                % For Pulseq version 1.4.0
                gx(j,i) = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:,i,j), lims);
                gy(j,i) = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:,i,j), lims);

                % ADC 
                adc = mr.makeAdc(adcSamples,lims,'Dwell',adcDwell);
%                 if params.spi.type == 0 
%                     adc = mr.makeAdc(adcSamples,lims,'Dwell',adcDwell);
%                 elseif params.spi.type == 1
%                     adc = mr.makeAdc(adcSamples,lims,'Dwell',adcDwell);
% %                     gx(j,i).delay = mr.calcDuration(gz_blips(end));
% %                     gy(j,i).delay = mr.calcDuration(gz_blips(end));
%                 end
                adc_post = mr.makeAdc(100,lims,'Dwell',adcDwell);
                % adc = mr.makeAdc(adcSamples,'Dwell',adcDwell);
%                 % Pulseq (and siemens) define the samples to happen in the center of the dwell period
%                 time_to_center = adc.dwell*((adcSamples-1)/2+0.5);
%                 adc.delay = round((gx.riseTime+gx.flatTime/2-time_to_center)/lims.rfRasterTime)*lims.rfRasterTime;

                % Getting the correct number to split the ADC
                % The ADC obj has to be splitted into N equal parts, with duration multiple of 10us
                adc_total_samples = adcSamples;
                % if params.gen.dork; adc_total_samples = adc_total_samples + adc_post.numSamples; end
                for k = 1:50
                    if mod(adc_total_samples,k) == 0 && mod(adc_total_samples/k*adcDwell,10e-9) == 0 && (adc_total_samples/k) < 8192 && mod(adc_total_samples/k,4) == 0
                        adcSplit = adc_total_samples/k;
                        break
                    end
                end
                params.gen.adc_split = adcSplit;
                params.gen.adc_segments = k;

                if params.spi.type == 0
%                     Readout ramp down
%                     AMM: need to nicely define the time (0.001) of this ramp gradients
                    gx_ramp(j,i) = mr.makeExtendedTrapezoid('x','times',[0 0.001],'amplitudes',[spiral_grad_shape(1,end,i,j),0]);
                    gy_ramp(j,i) = mr.makeExtendedTrapezoid('y','times',[0 0.001],'amplitudes',[spiral_grad_shape(2,end,i,j),0]);
                    if params.gen.echos > 1
                        tmp = cumsum(gx(j,i).waveform);
                        tmp = tmp./lims.gamma*100;   
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*6.5/100);
%                         gx_pre(i) = mr.makeTrapezoid('x',lims,'maxGrad',25e-3*lims.gamma,'Area',area);
                        gx_pre(i) = mr.makeTrapezoid('x',lims,'Area',area);
                        tmp = cumsum(gy(j,i).waveform);
                        tmp = tmp./lims.gamma*100;   
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*6.5/100);
%                         gy_pre(i) = mr.makeTrapezoid('y',lims,'maxGrad',25e-3*lims.gamma,'Area',area);
                        gy_pre(i) = mr.makeTrapezoid('y',lims,'Area',area);
                    end
                elseif params.spi.type == 1 || params.spi.type == 3
                    gx_ramp(j,i) = mr.makeExtendedTrapezoid('x','times',[0 0.0001],'amplitudes',[0,spiral_grad_shape(1,1,i,j)]);
                    gy_ramp(j,i) = mr.makeExtendedTrapezoid('y','times',[0 0.0001],'amplitudes',[0,spiral_grad_shape(2,1,i,j)]);

                    % Pre-phasing gradients, I adjust for 6.5% trajectory
                    % error to be closer to zero at center of k-space
                    % for rxy=3=6.5, rxy=4=
                    err = 6.5;
                    tmp = cumsum(gx(j,i).waveform);
                    tmp = tmp./lims.gamma*100;  
                    if params.spi.type == 1
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*err/100);
%                         area = ((tmp(1)-tmp(end))*4);
                    elseif params.spi.type == 3
                        area = ((tmp(1)-tmp(end/2))*4)+(((tmp(1)-tmp(end/2))*4)*err/100);
                    end
                    gx_pre = mr.makeTrapezoid('x',lims,'maxGrad',25e-3*lims.gamma,'Area',area);
%                     area = (max(abs(cumsum(gy(j,i) .waveform)))./lims.gamma*1000-abs((sum(gy_ramp(j,i).first:-1e4:gy_ramp(j,i).last)./1e3)))/2;
%                     area = (abs(min(cumsum(gy(j,i).waveform))) + abs(max(cumsum(gy(j,i).waveform))))/2./lims.gamma*1000;
                    tmp = cumsum(gy(j,i).waveform);
                    tmp = tmp./lims.gamma*100;   
                    if params.spi.type == 1
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*err/100);
%                         area = ((tmp(1)-tmp(end))*4);
                    elseif params.spi.type == 3
                        area = ((tmp(1)-tmp(end/2))*4)+(((tmp(1)-tmp(end/2))*4)*err/100);
                    end
                    gy_pre = mr.makeTrapezoid('y',lims,'maxGrad',25e-3*lims.gamma,'Area',area);

                end
            end
        end
elseif params.gen.ro_type == 'c'
    % Prepare EPI
    gy_pre_area = -(params.gen.del_k(2)*params.gen.n(2)/2/params.epi.ry)*params.epi.pf;
    if params.epi.pf ~= 1; gy_pre_area = gy_pre_area+(gy_pre_area/4*-1); end
    gy_pre = mr.makeTrapezoid('y',lims,'maxGrad',15e-3*lims.gamma,'Area',gy_pre_area);
    gy_blip = mr.makeTrapezoid('y',lims,'Area',params.gen.del_k(2));

    tot_bw = params.gen.n(1)*params.epi.bw_px;
    adcDwell = floor((1./tot_bw/lims.adcRasterTime))*lims.adcRasterTime;

    % ADC
    adcSamples = params.gen.n(1);
    % In SIEMENS number of ADC samples should be divisible by 4
    adcSamples = floor(adcSamples/4)*4;  

    % Making sure we are aligned to adcRastertime
    params.epi.bw = 1./adcDwell;
    adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',mr.calcDuration(gy_blip)/2); % Ramp sampling

    % Getting the correct number to split the ADC
    % The ADC obj has to be splitted into N equal parts, with duration multiple of 10us
    adc_total_samples = adcSamples*round(params.gen.n(2)/params.epi.ry*(params.epi.pf))*params.gen.n(3);
    % if params.gen.dork; adc_total_samples = adc_total_samples + adc_post.numSamples; end
    for i = 1:500
        if mod(adc_total_samples,i) == 0 && mod(adc_total_samples/i*adcDwell,10e-9) == 0 && (adc_total_samples/i) < 8192 && mod(adc_total_samples/i,4) == 0
            adcSplit = adc_total_samples/i;
            break
        end
    end
    params.gen.adc_split = adcSplit;
    params.gen.adc_segments = i;
    
    tmp = floor((adc.duration+(adc.delay*2))./lims.gradRasterTime)*lims.gradRasterTime;
    gx = mr.makeTrapezoid('x',lims,'Area',params.gen.n(1)*params.gen.del_k(1),'Duration',tmp);
    gx_pre = mr.makeTrapezoid('x',lims,'maxGrad',30e-3*lims.gamma,'Area',-gx.area/2);
    gx.amplitude = -gx.amplitude;

    % Pulseq (and siemens) define the samples to happen in the center of the dwell period
    time_to_center=adc.dwell*((adcSamples-1)/2+0.5);
    adc.delay = round(adc.delay/lims.rfRasterTime)*lims.rfRasterTime;

    % Let's split the gy_blip
    gy_blip_up = mr.makeExtendedTrapezoid('y',lims,'times',[0 gy_blip.riseTime],'amplitude',[0 gy_blip.amplitude]);
    gy_blip_down = mr.makeExtendedTrapezoid('y',lims,'times',[0 gy_blip.riseTime],'amplitude',[gy_blip.amplitude 0]);
    [gy_blip_up,gy_blip_down,~] = mr.align('right',gy_blip_up,'left',gy_blip_down,gx);
    gy_blips = mr.addGradients({gy_blip_down, gy_blip_up},lims);
end

% Flip Angles
% tr_tmp Rough estimate of TR, to calculate Ernst Angle
% n_lines = round(params.gen.n(2)/params.epi.ry)-((round(params.gen.n(2)/params.epi.ry)-round(params.gen.n(2)/params.epi.ry*params.epi.pf))/2);
n_lines = ceil(params.gen.n(2)/params.epi.ry)-((round(params.gen.n(2)/params.epi.ry)-round(params.gen.n(2)/params.epi.ry*params.epi.pf)));
if params.gen.ro_type == 'c'
    tr_tmp = (mr.calcDuration(gx)*(n_lines+3))+mr.calcDuration(gx_spoil)+params.gen.te; % +3 of navigators
elseif params.gen.ro_type == 's'
    tr_tmp = (mr.calcDuration(gx)*params.gen.echos)+2.6e-3+params.gen.te;
end

% Flip Angles
if params.gen.fa == 0
    params.gen.fa(1) = acos(exp(-(tr_tmp)/(params.gen.ernst_t1)));
else
    params.gen.fa(1) = params.gen.fa*pi/180;
end

% RF pulses depending on sequence type (VASO,ABC) and field strength,
% Values for: Duration, apodization and timeBwProduct, can be modified to
% avoid SAR or other issues ABC
if params.gen.seq == 2
    [rf0(1), gz(1)] = mr.makeSincPulse(params.gen.fa(1),'system',lims,'Duration',1e-3,...
        'SliceThickness',params.gen.fov(3),'apodization',0.5);
% VASO/BOLD/Fieldmap
else
    % 7T
    if params.gen.field_strength == 7 || params.gen.field_strength == 7i 
        [rf0(1), gz(1)] = mr.makeSincPulse(params.gen.fa(1),'system',lims,'Duration',2.56e-3,...
            'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
    % 9.4T
    elseif params.gen.field_strength == 9
        [rf0(1), gz(1)] = mr.makeSincPulse(params.gen.fa(1),'system',lims,'Duration',2.56e-3,...
            'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
%         [rf0(1), gz(1)] = mr.makeSincPulse(params.gen.fa(1),'system',lims,'Duration',5.25e-3,...
%             'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
    end
end
gzReph(1) = mr.makeTrapezoid('z',lims,'Area',-gz(1).area/2);


for i=2:params.gen.n(3)
    % Array of Flip Angles
    if params.gen.vfa
            e1 = exp(-tr_tmp/params.gen.ernst_t1);
            params.gen.fa(i) = asin((sin(params.gen.fa(1))*tan(params.gen.fa(i-1))) ...
                /((e1*sin(params.gen.fa(1)))+((1-e1)*tan(params.gen.fa(i-1)))));
    else
        if params.gen.fa == 0
            params.gen.fa(i) = acos(exp(-(tr_tmp)/(params.gen.ernst_t1)))*(180/pi);
        else
            params.gen.fa(i) = params.gen.fa(1);
        end
    end
    
    if ~isreal(params.gen.fa(i))
        params.gen.fa(i) = 90/180*pi;
    end
    
    % Array of RF pulses, depending on sequence type:
    % ABC
    if params.gen.seq == 2
        [rf0(i), gz(i)] = mr.makeSincPulse(params.gen.fa(i),'system',lims,'Duration',1e-3,...
            'SliceThickness',params.gen.fov(3),'apodization',0.5);
    % VASO/BOLD/Fieldmap
    else
        % 7T
        if params.gen.field_strength == 7 || params.gen.field_strength == 7i 
            [rf0(i), gz(i)] = mr.makeSincPulse(params.gen.fa(i),'system',lims,'Duration',2.56e-3,...
                'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
        % 9.4T
        elseif params.gen.field_strength == 9
            [rf0(i), gz(i)] = mr.makeSincPulse(params.gen.fa(i),'system',lims,'Duration',2.56e-3,...
                'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
        end
    end

    gzReph(i) = mr.makeTrapezoid('z',lims,'Area',-gz(i).area/2);
end

%% Prepare Delays and triggers
if params.gen.te > 0; te_delay = mr.makeDelay(round(params.gen.te-((mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph)+mr.calcDuration(gz_blips)),4)); end         % TE delay
if params.gen.tr_delay > 0; tr_delay = mr.makeDelay(params.gen.tr_delay); end                       % TR delay
if params.vaso.f_v_delay > 0; f_v_delay = mr.makeDelay(params.vaso.f_v_delay);          end         % FOCI-VASO delay
if params.vaso.v_b_delay > 0; v_b_delay = mr.makeDelay(params.vaso.v_b_delay);          end         % VASO-BOLD delay
if params.vaso.b_f_delay > 0; b_f_delay = mr.makeDelay(params.vaso.b_f_delay);          end         % BOLD-FOCI delay
% TR delay for Cartesian, to avoid PNS
if params.gen.ro_type == 'c'
    if params.epi.tr_delay > 0
        tr_delay_epi = mr.makeDelay(params.epi.tr_delay);  
    end
end
% Skope delays
if params.gen.skope ~= 0
    % AMM: Todo Need to confirm the times of this delays
    sk_pre_delay    = mr.makeDelay(1.84);
    sk_int_delay    = mr.makeDelay(200e-6);
    sk_post_delay   = mr.makeDelay(3);
    sk_min_tr_delay = gzReph(1).riseTime+gzReph(1).flatTime+gzReph(1).fallTime + ...
                      gz_blips(1).riseTime + gz_blips(1).flatTime + gz_blips(1).fallTime ...
                      + mr.calcDuration(gx(1)) + mr.calcDuration(gx_ramp(1));
    sk_min_tr_delay = mr.makeDelay(120e-3-sk_min_tr_delay);                         % Here I take Skope minTR = 110ms
    skope_trig = mr.makeDigitalOutputPulse('osc0','duration',10e-6,'system',lims);  % Skope trigger
end
% Fieldmap scan delays
if params.gen.seq == 3
    if params.epi.te > 0
        tmp = (mr.calcDuration(rf0(1))/2)+mr.calcDuration(gzReph(1))+mr.calcDuration(gz_blips(1))+mr.calcDuration(gx_pre);
        tmp = tmp + (round(round(params.gen.n(2)/params.epi.ry)/2)*(mr.calcDuration(gx(1)))+mr.calcDuration(gy_blip(1)));
        for i=1:length(params.epi.te)
            fm_te_delay(i) = mr.makeDelay(round(params.epi.te(i)-tmp,4));
        end
    end
end
% Delay for repetitions without MT pulse, only for ABC (seq=2)
if params.gen.seq == 2
    no_mt_delay = mr.makeDelay(mr.calcDuration(MT)+mr.calcDuration(gx_spoil));
end
% No FatSat delay...
% fs_delay = mr.makeDelay(mr.calcDuration(rf_fs)+mr.calcDuration(gx_fs));
% External trigger for fMRI
ext_trig = mr.makeDigitalOutputPulse('ext1','duration',10e-6,'system',lims);        % External trigger
dummy_delay = mr.makeDelay(10e-3);                                                  % Delay to use as dummy anywhere

%% Add blocks to Skope seq
if params.gen.skope == 1
    seq_sk=mr.Sequence();          % Create a new sequence object
    seq_sk.addBlock(sk_pre_delay); % Skope initial delay
    % Partitions loop
    for i=1:params.gen.n(3)
        % Spiral/Cartesian Interleaves/segments loop
        % for m=1:params.gen.n(3)
            for j = 1:params.gen.seg
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
                % AMM: skope here I only add the trigger in the center slice, (i=13)
                if i == 13
                    seq_sk.addBlock(skope_trig);
                    seq_sk.addBlock(sk_int_delay);     % Gradient free interval
                end
                % Gz blip
                if i < floor((params.gen.n(3)/2)+1)
    %                     seq.addBlock(gz_blips(i));
                        tmp_blip = gz_blips(i);
                elseif i > floor((params.gen.n(3)/2)+1)
    %                     seq.addBlock(gz_blips(i-1)); 
                        tmp_blip = gz_blips(i-1);
                end
%                 seq_sk.addBlock(skope_trig);
%                 seq_sk.addBlock(sk_int_delay);     % Gradient free interval
                if params.gen.ro_type == 's'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq_sk.addBlock(gx(i,j),gy(i,j),adc);
                    else
                        seq_sk.addBlock(gx(i,j),gy(i,j),tmp_blip,adc);
                    end
                    seq_sk.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                    % seq_sk.addBlock(gx_spoil,gy_spoil,gz_spoil);
                    if params.gen.tr_delay > 0; seq_sk.addBlock(tr_delay); end % TR delay
                elseif params.gen.ro_type == 'c'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq_sk.addBlock(gx_pre,gy_pre)
                    else
                        seq_sk.addBlock(gx_pre,gy_pre,tmp_blip)
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
                seq_sk.addBlock(sk_min_tr_delay); 
            end
        % end
    end
    seq_sk.addBlock(sk_post_delay);    % Skope post delay
end

%% Measure trajectory scan..
if params.gen.traj_scan
    seq_traj = mr.Sequence();
    
    % AMM: Need to define:
    % - slice thickness, 3mm or same as acq (1mm), - FA, -Duratin
    % 7T
    if params.gen.field_strength == 7
        [rf0_traj, gz_traj, gz_reph_traj] = mr.makeSincPulse(params.gen.fa(1),'system',lims,'Duration',3.5e-3,...
            'SliceThickness',params.gen.res(1)*0.75,'apodization',0.5,'timeBwProduct',25);
    % 9.4T
    elseif params.gen.field_strength == 9
        [rf0_traj, gz_traj, gz_reph_traj] = mr.makeSincPulse(params.gen.fa(1),'system',lims,'Duration',8.25e-3,...
            'SliceThickness',params.gen.res(1)*0.75,'apodization',0.5,'timeBwProduct',20);
    end
    
    % A few seconds delay...
    traj_scan_delay = mr.makeDelay(2000e-3);
    gx_traj = gz_traj;  gx_traj.channel = 'x';
    gy_traj = gz_traj;  gy_traj.channel = 'y';
    gx_reph_traj = gz_reph_traj(1);  gx_reph_traj.channel = 'x';
    gy_reph_traj = gz_reph_traj(1);  gy_reph_traj.channel = 'y';

    for i=1:params.spi.interl
        % Gx+ with readout gradient
        seq_traj.addBlock(rf0_traj,gx_traj);
        seq_traj.addBlock(gx_reph_traj);
        seq_traj.addBlock(gx(params.gen.n(3)/2+1,i),adc);
        seq_traj.addBlock(gx_ramp(1,1));
        seq_traj.addBlock(gz_spoil);
        seq_traj.addBlock(traj_scan_delay);
    
        % Gx- with readout gradient
        gx(params.gen.n(3)/2+1,i).waveform = gx(params.gen.n(3)/2+1,i).waveform.*(-1);
        seq_traj.addBlock(rf0_traj,gx_traj);
        seq_traj.addBlock(gx_reph_traj);
        seq_traj.addBlock(gx(params.gen.n(3)/2+1,i),adc);
        seq_traj.addBlock(gx_ramp(1,1));
        seq_traj.addBlock(gz_spoil);
        seq_traj.addBlock(traj_scan_delay);
    
        % Gx w/o readout gradient
        seq_traj.addBlock(rf0_traj,gx_traj);
        seq_traj.addBlock(gx_reph_traj);
        seq_traj.addBlock(adc);
        seq_traj.addBlock(gx_ramp(1,1));
        seq_traj.addBlock(gz_spoil);
        seq_traj.addBlock(traj_scan_delay);

        % Let's reset the sign of Gx
        gx(params.gen.n(3)/2+1,i).waveform  = gx(params.gen.n(3)/2+1,i).waveform.*(-1);
    end

    for i=1:params.spi.interl
        % Gy+ with readout gradient
        seq_traj.addBlock(rf0_traj,gy_traj);
        seq_traj.addBlock(gy_reph_traj);
        seq_traj.addBlock(gy(params.gen.n(3)/2+1,i),adc);
        seq_traj.addBlock(gy_ramp(1,1));
        seq_traj.addBlock(gz_spoil);
        seq_traj.addBlock(traj_scan_delay);
        
        % Gy- with readout gradient
        gy(params.gen.n(3)/2+1,i).waveform  = gy(params.gen.n(3)/2+1,i).waveform.*(-1);
        seq_traj.addBlock(rf0_traj,gy_traj);
        seq_traj.addBlock(gy_reph_traj);
        seq_traj.addBlock(gy(params.gen.n(3)/2+1,i),adc);
        seq_traj.addBlock(gy_ramp(1,1));
        seq_traj.addBlock(gz_spoil);
        seq_traj.addBlock(traj_scan_delay);
    
        % Gy w/o readout gradient
        seq_traj.addBlock(rf0_traj,gy_traj);
        seq_traj.addBlock(gy_reph_traj);
        seq_traj.addBlock(adc);
        seq_traj.addBlock(gy_ramp(1,1));
        seq_traj.addBlock(gz_spoil);
        seq_traj.addBlock(traj_scan_delay);

        % Let's reset the sign of Gy
        gy(params.gen.n(3)/2+1,i).waveform  = gy(params.gen.n(3)/2+1,i).waveform.*(-1);
    end

end

%% Actual scan
seq = mr.Sequence();
seq.addBlock(ext_trig);                                      % External trigger
if params.gen.seq == 1
    if params.vaso.foci; seq.addBlock(rf_foci); end              % FOCI
    if params.vaso.f_v_delay > 0; seq.addBlock(f_v_delay); end   % FOCI-VASO delay
end

% if params.gen.seq == 2
%     if params.mt.mt == 1
%         % if params.gen.fat_sat; seq.addBlock(rf_fs,gx_fs,gy_fs,gz_fs);   end      % fat-sat
%         seq.addBlock(MT);                                           % MT pulse 
%         seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
%     end
% end

% dur0 = 0;
% dur1 = 0;
last_part = params.gen.n(3);
% if params.gen.echos > 1 && params.spi.interl > 1; last_part = last_part/2; end
for i_ro_blocks = 1:ro_blocks
    for i=1:last_part
%         if params.gen.fat_sat; seq.addBlock(rf_fs,gx_fs,gy_fs,gz_fs);   end      % fat-sat
        % Adding MT pulse every specified time
        if params.gen.seq == 2
            if params.mt.mt == 1 
                    seq.addBlock(MT);                                      % MT pulse 
                    seq.addBlock(gz_spoil);

            else
                    seq.addBlock(no_mt_delay);
            end
        end

        for j=1:params.gen.seg
            if params.gen.fat_sat; seq.addBlock(rf_fs,gx_fs,gy_fs,gz_fs);   end      % fat-sat
%             if params.gen.fat_sat && j>1; seq.addBlock(fs_delay);   end        % fat-sat delayy...
            if i==1 && j==1 && i_ro_blocks==1; tr0 = seq.duration(); end                           % save seq dur to calc TR
            rf = rf0(i);
            rf.phaseOffset = rf_phase_offset(i);
            adc.phaseOffset = adc_phase_offset(i);
            adc_post.phaseOffset = adc_phase_offset(i);
            seq.addBlock(rf,gz(i));
            if i==1; te0 = seq.duration(); end                           % save seq dur to calc TE
            seq.addBlock(gzReph(i));
            % EPI navigators
            if params.gen.ro_type == 'c'
                for i_nav = 1:3
                    gx.amplitude = -gx.amplitude;
                    seq.addBlock(gx,adc); 
                end
            end
            % Gz blip
            if i < floor((params.gen.n(3)/2)+1) && params.gen.kz_enc == 0
%                     seq.addBlock(gz_blips(i));
                    tmp_blip = gz_blips(i);
            elseif or(i > floor((params.gen.n(3)/2)+1) && params.gen.kz_enc == 0, i > 1 && params.gen.kz_enc == 1)
%                     seq.addBlock(gz_blips(i-1)); 
                    tmp_blip = gz_blips(i-1);
            end
            % Echos loop
            for k=1:params.gen.echos
                if params.gen.ro_type == 's'
                    if k==1 || k==params.gen.echos-1
                        if or(i == floor((params.gen.n(3)/2)+1) && params.gen.kz_enc == 0, i == 1 && params.gen.kz_enc == 1)
                            if params.gen.skope == 2
                                seq.addBlock(skope_trig);
                                seq.addBlock(sk_int_delay);     % Gradient free interval
                            end
                        else
                            seq.addBlock(tmp_blip);
                        end
                    end
                    if and(params.gen.seq == 3, params.epi.te > 0); seq.addBlock(fm_te_delay(i_ro_blocks)); end
                    if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
                    % if spiral-in or in-out, we need in-plane pre-phasing
                    if params.spi.type == 1 || params.spi.type == 3 
                        if params.spi.interl > 1
                            gx_pre.amplitude = gx_pre.amplitude*-1;
                            gy_pre.amplitude = gy_pre.amplitude*-1;
                        end
                        seq.addBlock(gx_pre,gy_pre);                       
                    end
                    if i==1; te1 = seq.duration(); end                           % save seq dur to calc TE
                    seq.addBlock(gx(i,j),gy(i,j),adc);
                    seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                    if params.gen.dork; seq.addBlock(adc_post); end
                    if k==params.gen.echos; seq.addBlock(gz_spoil); end
                    % rewinder if multiple echos
                    if params.gen.echos > 1; seq.addBlock(gx_pre(j),gy_pre(j)); end
                    if params.gen.tr_delay > 0; seq.addBlock(tr_delay); end % TR delay
                elseif params.gen.ro_type == 'c'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq.addBlock(gx_pre,gy_pre)
                    else
                        seq.addBlock(gx_pre,gy_pre,tmp_blip)
                    end
                    for k = 1:n_lines
                        gx.amplitude = -gx.amplitude;
                        if k == 1
                            seq.addBlock(gx,gy_blip_up,adc);
                        elseif k == n_lines
                            seq.addBlock(gx,gy_blip_down,adc);
                        else
                            seq.addBlock(gx,gy_blips,adc);
                        end
                        if i==1 && k==ceil((n_lines/2)+1); te1 = seq.duration(); end                           % save seq dur to calc TE
                    end
                    if gx.amplitude > 0
                        gx.amplitude = -gx.amplitude;
                    end
                    seq.addBlock(gz_spoil);
                    if params.gen.tr_delay > 0; seq.addBlock(tr_delay); end % TR delay
                    if params.epi.tr_delay > 0; seq.addBlock(tr_delay_epi); end
                end
            end
            if i==1 && j==1 && i_ro_blocks==1; tr1 = seq.duration(); end                           % save seq dur to calc TR
%             % Adding MT pulse every specified time
%             if params.gen.seq == 2
%                 dur1 = seq.duration();
%                 if params.mt.mt == 1 && i_ro_blocks == 1
%                     if dur1-dur0 > params.mt.mt_rep  && j == params.gen.seg
% %                         if params.gen.fat_sat; seq.addBlock(rf_fs,gx_fs,gy_fs,gz_fs);   end      % fat-sat
%                         seq.addBlock(MT);                                      % MT pulse 
%                         seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
%                         dur0 = dur1;
%                     else
% %                         if params.gen.fat_sat; seq.addBlock(rf_fs,gx_fs,gy_fs,gz_fs);   end      % fat-sat
%                         seq.addBlock(no_mt_delay);
%                     end
%                 else
% %                         if params.gen.fat_sat; seq.addBlock(rf_fs,gx_fs,gy_fs,gz_fs);   end      % fat-sat
%                         seq.addBlock(no_mt_delay);
%                 end
%             end

    %         seq.addBlock(dummy_delay);          % AMM: Temp: Adding a delay to let mag recover
        end
        dur1 = seq.duration();      
    end
    if params.gen.seq == 1 && i_ro_blocks == 1
        if params.vaso.v_b_delay > 0; seq.addBlock(v_b_delay); end      % VASO-BOLD delay
    elseif params.gen.seq == 1 &&  i_ro_blocks == 2
        if params.vaso.b_f_delay > 0; seq.addBlock(b_f_delay); end   % BOLD-FOCI delay
    end
end

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
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
j = 1;
if params.gen.ro_type == 's'
    plane_samples = adc.numSamples;
elseif params.gen.ro_type == 'c'
    plane_samples = adc.numSamples*n_lines;
    % Discarding the EPI navigator samples, here I have 3
    j = j+(adc.numSamples*3);
    tmp = ktraj_adc(1,j:j+adc.numSamples*3);
    tmp_mx = max(tmp(:));
end
for i=1:params.gen.n(3)
    l = 1;
    for m=1:params.spi.interl
            ks_traj.kx(l:l+plane_samples-1,i) = ktraj_adc(1,j:j+plane_samples-1);
            ks_traj.ky(l:l+plane_samples-1,i) = ktraj_adc(2,j:j+plane_samples-1);
            ks_traj.kz(l:l+plane_samples-1,i) = ktraj_adc(3,j:j+plane_samples-1);
            j = j+(plane_samples*params.gen.echos);
            l = l+plane_samples;
            if params.gen.dork; j = j+adc_post.numSamples;  end
            if params.gen.ro_type == 'c'
                j = j+(adc.numSamples*3);
            end
    end
end
% Splitting trajectories for ech1 and ech2 of IN-OUT
% AMM ToDo : Fix this part so it works for more echos..
if params.spi.type == 3
    % Echo 1 (IN) trajectory
    ks_traj_e1.kx = ks_traj.kx([1:(plane_samples/2)-1,plane_samples+1:end-(plane_samples/2)-1],:);
    ks_traj_e1.ky = ks_traj.ky([1:(plane_samples/2)-1,plane_samples+1:end-(plane_samples/2)-1],:);
    ks_traj_e1.kz = ks_traj.kz([1:(plane_samples/2)-1,plane_samples+1:end-(plane_samples/2)-1],:);
    ks_traj_e1.kx = [padarray(ks_traj_e1.kx(1:end/2,:),[1 0],'post');padarray(ks_traj_e1.kx((end/2)+1:end,:),[1 0],'post')];
    ks_traj_e1.ky = [padarray(ks_traj_e1.ky(1:end/2,:),[1 0],'post');padarray(ks_traj_e1.ky((end/2)+1:end,:),[1 0],'post')];
    ks_traj_e1.kz = [padarray(ks_traj_e1.kz(1:end/2,:),[1 0],'post');padarray(ks_traj_e1.kz((end/2)+1:end,:),[1 0],'post')];
     % Echo 2 (OUT) trajectory
    ks_traj_e2.kx = ks_traj.kx([(plane_samples/2):plane_samples-2,plane_samples+(plane_samples/2):end-2],:);
    ks_traj_e2.ky = ks_traj.ky([(plane_samples/2):plane_samples-2,plane_samples+(plane_samples/2):end-2],:);
    ks_traj_e2.kz = ks_traj.kz([(plane_samples/2):plane_samples-2,plane_samples+(plane_samples/2):end-2],:);
    ks_traj_e2.kx = [padarray(ks_traj_e2.kx(1:end/2,:),[1 0],'pre');padarray(ks_traj_e2.kx((end/2)+1:end,:),[1 0],'pre')];
    ks_traj_e2.ky = [padarray(ks_traj_e2.ky(1:end/2,:),[1 0],'pre');padarray(ks_traj_e2.ky((end/2)+1:end,:),[1 0],'pre')];
    ks_traj_e2.kz = [padarray(ks_traj_e2.kz(1:end/2,:),[1 0],'pre');padarray(ks_traj_e2.kz((end/2)+1:end,:),[1 0],'pre')];
end

% % Plotting traj partition by partition
figure(17);
plot3(ks_traj.kx(:),ks_traj.ky(:),ks_traj.kz(:),'DisplayName',sprintf('Kz = %i',i)); title('3D K-space'); xlabel('Kx (1/m)'); ylabel('Ky (1/m)'); zlabel('Kz (1/m)');
hold on
view(2)

% % Plotting trajectory by segment
% figure(15);
% hold on
% for i = 1:params.gen.n(3)
%     idx=1;
%     for j=1:params.spi.interl                                    
%         plot3(ks_traj.kx(idx:idx+adcSamples-1,i),ks_traj.ky(idx:idx+adcSamples-1,i),ks_traj.kz(idx:idx+adcSamples-1,i)); title('3D K-space');
%         idx = idx+adcSamples;
%     end
% end
% view(2)

%% Adding some extra parameters to params directory
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
params.gen.n(1:2) = params.gen.n(1:2)+tmp; % making in-plane even
% Adjust for phase oversampling
params.gen.n_ov(1:2) = params.gen.n(1:2);
params.gen.n_ov(3) = params.gen.n_ov(3)/params.gen.kz;
[params.gen.n, params.gen.n_ov] = deal(params.gen.n_ov, params.gen.n);

% converting FA back to degrees
params.gen.fa = real(params.gen.fa)*180/pi;

%% Adding ro_samples, acqTR,volTR, effTR, TE and other params
params.gen.ro_samples = adc.numSamples;
params.gen.TE = te1-te0+(mr.calcDuration(rf)/2);
params.gen.acqTR = tr1-tr0;
params.gen.volTR = params.gen.acqTR.*params.gen.n(3).*params.spi.interl;
params.gen.effTR = seq.duration();
params.gen.adc_dwell = adc.dwell;
params.gen.ro_time = adc.duration;
params.gen.ti1 = ((params.gen.effTR-params.vaso.f_v_delay)/4)+params.vaso.f_v_delay;
if params.spi.type ==3; params.gen.echos = 2; end       % If IN-OUT.. 2 echos

%% Calculating time vector for MRIReco.jl
% AMM: Here I might need to take into account the gradient delays...
% julia_time = params.gen.TE+adc.delay+(adc.dwell:adc.dwell:adc.duration);
if params.gen.ro_type == 's'
    julia_time = repmat(params.gen.TE+(0:adc.dwell:adc.duration-adc.dwell),1,params.spi.interl);
elseif params.gen.ro_type == 'c'
    julia_time = params.gen.TE+(0:adc.dwell:(adc.duration*n_lines)-adc.dwell);
end
params.gen.t_vector = julia_time;

% Adjusting for different readout types...
if params.spi.type == 1
    params.gen.TE = params.gen.TE+mr.calcDuration(gx);
elseif params.spi.type == 3
    params.gen.TE = params.gen.TE+(mr.calcDuration(gx)/2);
    params.gen.ro_samples = params.gen.ro_samples/2;
    params.gen.t_vector = params.gen.t_vector([1:plane_samples/2,plane_samples+1:(plane_samples)+(plane_samples/2)]);
end

%% Check accoustic resonance frequencies, Des script
% t     time of samples in microseconds, y     sample values, Fs    sampling rate in Hz
Fs = 1/lims.adcRasterTime/100;
% gradients = seq.gradient_waveforms1()/lims.gamma*1000;   % Full sequence
% gradients = gradients(1:2,:);               % taking only gx and gy
% time = linspace(0,seq.duration(),size(gradients,2)).*1e6; % Full sequence
gradients = [gx(1).waveform; gy(1).waveform]./lims.gamma*1000;   % 1 readout
time = linspace(0,mr.calcDuration(gx(1)),size(gradients,2)).*1e6; % 1 readout
gaxes = ['X' 'Y'];
for i=1:length(gaxes)   
    gradFreqPlot_pulseq(time,gradients(i,:),Fs,gaxes(i),params.gen.field_strength);
end

%% Set definitions
seq.setDefinition('MaxAdcSegmentLength',params.gen.adc_split);
seq.setDefinition('FOV', params.gen.fov);
seq.setDefinition('Name', seq_name);

%% Saving files
% Check if folder exist, if no it creates it
if exist(sprintf('./data/%s',folder_name)) == 0
    tmp = sprintf('./data/%s',folder_name);
    system(sprintf('mkdir %s',tmp));
    system(sprintf('mkdir %s/acq',tmp));
%     system(sprintf('mkdir %s/acq/romeo',tmp));
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
if params.gen.skope == 1
    seq_sk.write(strcat(namestr,'_sk.seq'));
end
if params.gen.traj_scan == 1
    seq_traj.write(strcat(namestr,'_traj_scan.seq'));
end
if params.spi.type == 3
    % AMM ToDo : Fix this part so it works for more echos..
    ks_traj = ks_traj_e1;
    save(strcat(namestr,'_e1_ks_traj_nom.mat'),'ks_traj')
    ks_traj = ks_traj_e2;
    save(strcat(namestr,'_e2_ks_traj_nom.mat'),'ks_traj')
else
    save(strcat(namestr,'_ks_traj_nom.mat'),'ks_traj')
end
save(strcat(namestr,'_params.mat'),'params')
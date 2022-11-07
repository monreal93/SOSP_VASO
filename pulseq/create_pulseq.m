clear all; clc
% Change directory to path where sosp_vaso git has been cloned
cd /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/
addpath(genpath("./pulseq/functions"))
% Add path to pulseq installation folder
addpath(genpath("/home/amonreal/Documents/PhD/tools/pulseq/"))
addpath(genpath("/home/amonreal/Documents/PhD/tools/tOptGrad_V0.2/minTimeGradient/mex_interface/"))
addpath(genpath("/home/amonreal/Documents/PhD/PhD_2022/BSc_students/"))

%% ToDo
% - Check exactly how the saturation pulse should be, what spoliers I need?
% - How am I gonna take into account the gradient delays (trajectory,adc
% and time vector)
% - ABC: implement partition encoding direction center-out and then segment
% interleaved (Viktor's email)
% - check values of duration, timeBwProduct rf0
% - make params.gen.pf work...
% - Include epi.pf in volTR estimation for Cartesian
% - Include Gauss pulse in volTR for ABC seq
% - Fix TE calculation for Cartesian
% - Allow segmentation of EPI train...
% - Check that volTR takes into acount kz
% - Check update of n/res with 'c' and 'pf' seems wrong, slices gets 22
% - Properly define adcDwell time for gx in 'c' (now I use 2.2ms) 'FlatTime'
% - check values of duration, timeBwProduct, rf0
% - Properly define golden angle in prepare_spirals_rf_grad_adc.m function
% - Define if in tr_tmp, I want to include the time from fat sat..for FA

%% Define parameters
% Set system limits (MaxGrad=70, MaxSlew = 200);
% lims = mr.opts('MaxGrad',65,'GradUnit','mT/m',...
%     'MaxSlew',190,'SlewUnit','T/m/s',...
%     'rfRingdownTime', 30e-6,'rfDeadtime', 180e-6,'adcDeadTime', 10e-6, 'B0',7);  % To read it in VM I need rfDeadtime = 180e-6
lims = mr.opts('MaxGrad',65,'GradUnit','mT/m',...
    'MaxSlew',190,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6,'rfDeadtime', 180e-6,'adcDeadTime', 10e-6, 'B0',6.98);  % To read it in VM I need rfDeadtime = 180e-6

folder_name = '11012022_abc';         % Day I am scanning
seq_name = 'sample';                        % use sv_n (n for the diff scans at each day)
params.gen.seq = 2;                         % 1-VASO 2-ABC 3-Fieldmap

% General parameters
params.gen.fov = [192 192 24].*1e-3;
params.gen.res = [1.4 1.4 1].*1e-3;     % Target resoultion (for 0.8 use 0.78)
params.gen.fa = 0;                 % Set to 0, to use Erns Angle
params.gen.ernst_t1 = 1220e-3;      % T1 to calculate Ernst Angle (s)
params.gen.te = 0e-3;               % Set to 0 to shortest TE possible
params.gen.ro_type = 's';           % 's'-Spiral, 'c'-Cartesian
params.gen.kz = 1;                  % Acceleration in Kz
params.gen.pf = 1;                  % Partial fourier in Kz
params.gen.fat_sat = 0;             % Fat saturation (1=yes,0=no)
params.gen.skope = 0;               % Add skope sync scan and triggers
         
% Spiral parameters
params.spi.type = 1;                % spiral type 0=spiral-Out , 1=spiral-In
params.spi.rotate = 'none';         % Spiral rotation ('none','linear','golden'), linear not implemented
params.spi.increment = 'linear';    % Spiral increment mode (for now only linear)
params.spi.max_grad  = 35;          % Peak gradient amplitude for spiral (mT/m)  (35)
params.spi.max_sr = 155;            % Max gradient slew rate for spiral (mT/m/ms) (155)
params.spi.interl = 1;              % Spiral interleaves
params.spi.vd = 1.6;                % Variability density
params.spi.rxy = 3;                 % In-plane undersampling
params.spi.bw = 500e3;              % Spiral BW in MHz (Max value 1,000e3) (500e3)

% MT pulse parameters
params.mt.mt = 0;                   % Add MT pulse, 0 for reference scan without MT
params.mt.alpha = 225;  %225
params.mt.delta = 650;  %650      % for Pulseq approach should be 650 to match Viktors phase
params.mt.trf = 0.004;
params.mt.mt_rep = 130e-3;           % How often to play the mt pulse

% EPI parameters
params.epi.ry = 3;
params.epi.pf = 6/8;
params.epi.te = [33.6 36.6 38.6 40]*1e-3+0.022; % Echo times for field map
params.epi.seg = 1;                 % EPI Segments
params.epi.tr_delay = 20e-3;        % Delay after each TR, needed to reduce SAR
params.epi.bw_px = 960;             % BW in Hz  960/396

% VASO parameters
params.vaso.foci = 1;               % FOCI inversion?
params.vaso.tr = 4500e-3;           % volume TR
params.vaso.ti1 = 1800e-3;          % VASO TI1, 
params.vaso.ti2 = params.vaso.ti1+(params.vaso.tr/2);
params.vaso.foci_ampl = 270;        % FOCI amplitude (140/270)
params.vaso.foci_dur = 10410e-6;    % FOCI duration
params.vaso.foci_bw = 150;          % FOCI Bandwidth (150)
params.vaso.f_v_delay = 900e-3;     % FOCI-VASO delay
params.vaso.v_b_delay = 10e-3;      % VASO-BOLD delay
params.vaso.b_f_delay = 5e-3;       % BOLD-FOCI delay

% Some calculations, restrictions in seq...
params.gen.del_k = (1./params.gen.fov)*params.gen.kz;
params.gen.n = round(params.gen.fov./params.gen.res);

% Trying to make n multiple of 4,update res
tmp = mod(params.gen.n,4);
params.gen.n(1:2) = params.gen.n(1:2)+tmp(1:2);
% ToDo: Check... Not sure if I do need this...
% params.gen.n(3) = params.gen.n(3)/params.gen.kz/params.gen.pf;
if params.gen.seq == 3; params.gen.ro_type = 'c'; end   % if fieldmap, cartesian
if params.gen.seq == 2; params.gen.ro_type = 's'; end   % if ABC, Spiral
if params.gen.ro_type == 'c'
    params.gen.seg = params.epi.seg;
    % Trying to make n a multiple of 
elseif params.gen.ro_type == 's'
    params.gen.seg = params.spi.interl;
end

%% Create blocks (I will get this into functions)
% TR-FOCI
[B1_foci,phase,rf_complex,Gz_foci,fa] = tr_foci(params);
% Emily's implementation:
% HSIR.m_dFlipAngle/HSIR.m_sFociData.FA_scale  240/0.6512
% fa = 240/0.6512;
rf_foci = mr.makeArbitraryRf(rf_complex,fa*pi/180, 'system', lims);%, 'Delay',2e-4);

% MT pulse
if params.gen.seq == 2
%     % Using Viktor's code
%     vpulse = VPF_gaussian_pulse_4_Maastricht(params.mt.alpha,params.mt.delta,params.mt.trf);
%     MT = mr.makeArbitraryRf(vpulse.b1,params.mt.alpha*pi/180, 'system', lims);%, 'Delay',2e-4);
    % Using Pulseq
    MT = mr.makeGaussPulse(params.mt.alpha*pi/180,lims,'Duration',params.mt.trf,'FreqOffset',params.mt.delta);
end

% Fat sat:
% AMM: ToDo: how to properly desing the fat-sat pulse
sat_ppm = -3.45;
fs_angle = 110;
sat_freq = sat_ppm*1e-6*lims.B0*lims.gamma;
rf_fs = mr.makeGaussPulse(fs_angle*pi/180,'system',lims,'Duration',5e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
gx_fs = mr.makeTrapezoid('x',lims,'delay',mr.calcDuration(rf_fs),'Area',-1/1e-4); % spoil up to 0.1mm
gy_fs = mr.makeTrapezoid('y',lims,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
gz_fs = mr.makeTrapezoid('z',lims,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

%% Prepare kz Blips:
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
        gz_blips(i) = mr.makeTrapezoid('z',lims,'Area',area);
    end
end
if mod(floor(params.gen.n(3)/params.gen.kz),2) == 1
    gz_blips(round(params.gen.n(3)/2)) = [];    % Removing the empty blip...
else
    gz_blips(round(params.gen.n(3)/2)+1) = [];    % Removing the empty blip...
end

%% Preparing readout elements
[rf_phase_offset,adc_phase_offset] = rf_adc_phase(params);
if params.gen.ro_type == 's'
    [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc(params,lims);
        for j=1:params.gen.n(3)
            for i=1:params.spi.interl
                % Readout gradients
                gx(j,i) = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:,i,j), lims);
                gy(j,i) = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:,i,j), lims);
        
                % ADC
                % adcSamples = adcSamples - mr.calcDuration(gz_blips(1))/adcDwell;
%                 adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',mr.calcDuration(gz_blips(1)));
                adc = mr.makeAdc(adcSamples,'Dwell',adcDwell);
%                 % Pulseq (and siemens) define the samples to happen in the center of the dwell period
%                 time_to_center = adc.dwell*((adcSamples-1)/2+0.5);
%                 adc.delay = round((gx.riseTime+gx.flatTime/2-time_to_center)/lims.rfRasterTime)*lims.rfRasterTime;
                if params.spi.type == 0
%                     Readout ramp down
%                     AMM: need to nicely define the time (0.001) of this ramp gradients
                    gx_ramp(j,i) = mr.makeExtendedTrapezoid('x','times',[0 0.001],'amplitudes',[spiral_grad_shape(1,end,i,j),0]);
                    gy_ramp(j,i) = mr.makeExtendedTrapezoid('y','times',[0 0.001],'amplitudes',[spiral_grad_shape(2,end,i,j),0]);
                elseif params.spi.type == 1
                    gx_ramp(j,i) = mr.makeExtendedTrapezoid('x','times',[0 0.001],'amplitudes',[0,spiral_grad_shape(1,1,i,j)]);
                    gy_ramp(j,i) = mr.makeExtendedTrapezoid('y','times',[0 0.001],'amplitudes',[0,spiral_grad_shape(2,1,i,j)]);
                    area = (params.gen.del_k(1)*(params.gen.n(1)/2));
                    gx_pre = mr.makeTrapezoid('x',lims,'maxGrad',25e-3*lims.gamma,'Area',area);
                    area = (params.gen.del_k(2)*(params.gen.n(2)/2));
                    gy_pre = mr.makeTrapezoid('y',lims,'maxGrad',25e-3*lims.gamma,'Area',area*2);
                end
            end
        end
elseif params.gen.ro_type == 'c'
    % Prepare EPI

    gy_pre = mr.makeTrapezoid('y',lims,'maxGrad',20e-3*lims.gamma,'Area',-(params.gen.del_k(2)*params.gen.n(2)/2)*(params.epi.pf)...
        +(params.gen.del_k(2)*params.gen.n(2)/2)-(params.gen.del_k(2)*params.gen.n(2)/2)*(params.epi.pf));    
    gy_blip = mr.makeTrapezoid('y',lims,'Area',params.gen.del_k(2).*params.epi.ry);

    % AMM: Need to properly define the dwell time, now I have 2e-6
    % adcDwell = 2e-6;   % AMM: Here i use 2e-6 need to properly define it
    tot_bw = params.gen.n(1)*params.epi.bw_px;
    adcDwell = floor((1./tot_bw/lims.adcRasterTime))*lims.adcRasterTime;

    % ADC
    adcSamples = params.gen.n(1);
    % In SIEMENS number of ADC samples should be divisible by 4
    adcSamples = floor(adcSamples/4)*4;  

    % Making sure we are aligned to adcRastertime
%     tmp = (adcDwell*adcSamples)+(mr.calcDuration(gy_blip)/2);
    params.epi.bw = 1./adcDwell;

%     gx = mr.makeTrapezoid('x',lims,'riseTime',1e-4,'system',lims,'flatArea',params.gen.n(1)*params.gen.del_k(1),'FlatTime',params.gen.n(1)*2e-6); % Dwell time 2e-6
%     % No ramp sampling
%     % gx = mr.makeTrapezoid('x',lims,'flatArea',params.gen.n(1)*params.gen.del_k(1),'FlatTime',round(params.gen.n(1)*adcDwell,4)); % Dwell time 2.2e-6
    % Ramp sampling
%     gx = mr.makeTrapezoid('x',lims,'Area',params.gen.n(1)*params.gen.del_k(1));%,'Duration',round(params.gen.n(1)*adcDwell,4)); % Dwell time 2.2e-6
%     gx_pre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2);
%     gx.amplitude = -gx.amplitude;


%     adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',gx.riseTime-(tmp/2)); % No ramp sampling
    adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',mr.calcDuration(gy_blip)/2); % Ramp sampling

    % Getting the correct number to split the ADC
    % The ADC obj has to be splitted into N equal parts, with duration multiple
    % of 10us
    adc_total_samples = adcSamples*round(params.gen.n(2)/params.epi.ry*(params.epi.pf))*params.gen.n(3);
    for i = 1:100
        if mod(adc_total_samples,i) == 0 && mod(adc_total_samples/i*adcDwell,10e-9) == 0 && (adc_total_samples/i) < 8192 && mod(adc_total_samples/i,4) == 0
            adcSplit = adc_total_samples/i;
            break
        end
    end
    params.gen.adc_split = adcSplit;
    params.gen.adc_segments = i;
    
    tmp = floor((adc.duration+(adc.delay*2))./lims.gradRasterTime)*lims.gradRasterTime;
    gx = mr.makeTrapezoid('x',lims,'Area',params.gen.n(1)*params.gen.del_k(1),'Duration',tmp);%,'Duration',round(params.gen.n(1)*adcDwell,4)); % Dwell time 2.2e-6
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

% RF excitation
% ToDo: check values of duration, timeBwProduct
% Rough estimate of TR, to calculate Ernst Angle
if params.gen.ro_type == 'c'
    tr_tmp = mr.calcDuration(gx)*(round(params.gen.n(2)/params.epi.ry*(params.epi.pf))+3)+2.6e-3;
%     if params.gen.fat_sat == 1
%         tr_tmp = tr_tmp + mr.calcDuration(rf_fs) +  mr.calcDuration(gz_fs);
%     end
elseif params.gen.ro_type == 's'
    tr_tmp = mr.calcDuration(gx)+2.6e-3;
%     if params.gen.fat_sat == 1
%         tr_tmp = tr_tmp + mr.calcDuration(rf_fs) +  mr.calcDuration(gz_fs);
%     end
end

%% Flip Angles
if params.gen.fa == 0
    params.gen.fa = acos(exp(-(tr_tmp)/(params.gen.ernst_t1)))*(180/pi);
end

if params.gen.seq == 2
    [rf0, gz] = mr.makeSincPulse(params.gen.fa*pi/180,'system',lims,'Duration',1e-3,...
        'SliceThickness',params.gen.fov(3),'apodization',0.5);
else
    [rf0, gz] = mr.makeSincPulse(params.gen.fa*pi/180,'system',lims,'Duration',2.6e-3,...
        'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
end
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2);


%% Prepare spoilers
% AMM: Todo: Need to confirm the size (amplitude/area) of this spoiler
mag_spoil = 26e-3;     % mT
sr_spoil = 100;        % mT/m/s
gx_spoil=mr.makeTrapezoid('x','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);%,'delay',1e-3);
gy_spoil=mr.makeTrapezoid('y','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);%,'delay',1e-3);
gz_spoil=mr.makeTrapezoid('z','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);%,'delay',1e-3);

%% Prepare Delays, triggers
if params.gen.te > 0; te_delay = mr.makeDelay(round(params.gen.te-(max(rf0.t)-min(rf0.t))/2,4)); end         % TE delay
if params.vaso.f_v_delay > 0; f_v_delay = mr.makeDelay(params.vaso.f_v_delay);          end         % FOCI-VASO delay
if params.vaso.v_b_delay > 0; v_b_delay = mr.makeDelay(params.vaso.v_b_delay);          end         % VASO-BOLD delay
if params.vaso.b_f_delay > 0; b_f_delay = mr.makeDelay(params.vaso.b_f_delay);          end         % BOLD-FOCI delay
if params.gen.ro_type == 'c'; tr_delay = mr.makeDelay(params.epi.tr_delay);             end         % TR delay for Cartesian, to avoid PNS
if params.gen.skope == 1
    % AMM: Todo Need to confirm the times of this delays
    sk_pre_delay    = mr.makeDelay(1.84);
    sk_int_delay    = mr.makeDelay(200e-6);
    sk_post_delay   = mr.makeDelay(3);
    sk_min_tr_delay = gzReph.riseTime+gzReph.flatTime+gzReph.fallTime + ...
                      gz_blips(1).riseTime + gz_blips(1).flatTime + gz_blips(1).fallTime ...
                      + gx(1).shape_dur + gx_ramp(1).shape_dur;
    sk_min_tr_delay = mr.makeDelay(120e-3-sk_min_tr_delay);                         % Here I take Skope minTR = 110ms
    skope_trig = mr.makeDigitalOutputPulse('osc0','duration',10e-6,'system',lims);  % Skope trigger
end
% Delays for Fieldmap scan
if params.gen.seq == 3
    if params.epi.te > 0
        tmp = (mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph)+mr.calcDuration(gz_blips(1))+mr.calcDuration(gx_pre);
        tmp = tmp + (round(round(params.gen.n(2)/params.epi.ry)/2)*(mr.calcDuration(gx(1)))+mr.calcDuration(gy_blip(1)));
        for i=1:length(params.epi.te)
            fm_te_delay(i) = mr.makeDelay(round(params.epi.te(i)-tmp,4));
        end
    end
end
ext_trig = mr.makeDigitalOutputPulse('ext1','duration',10e-6,'system',lims);        % External trigger
dummy_delay = mr.makeDelay(10e-3);           % AMM: Temp: Delay to use as dummy anywhere

%% Add blocks to Skope seq
if params.gen.skope == 1
    seq_sk=mr.Sequence();          % Create a new sequence object
    seq_sk.addBlock(sk_pre_delay); % Skope initial delay
    % Partitions loop
    for i=1:params.gen.n(3)
        % Spiral/Cartesian Interleaves/segments loop
        for m=1:params.gen.n(3)
            for j = 1:params.gen.seg
                adc.phaseOffset = adc_phase_offset(i);
                seq_sk.addBlock(gzReph);
                if params.gen.te > 0; seq_sk.addBlock(te_delay); end % TE delay
                % EPI navigators
                if params.gen.ro_type == 'c'
                    for i_nav = 1:3
                        gx.amplitude = -gx.amplitude;
                        seq.addBlock(gx,adc); 
                    end
                end
                seq_sk.addBlock(skope_trig);
                seq_sk.addBlock(sk_int_delay);     % Gradient free interval
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
                        seq.addBlock(gx(i,j),gy(i,j),adc);
                    else
                        seq.addBlock(gx(i,j),gy(i,j),tmp_blip,adc);
                    end
                    seq_sk.addBlock(gx_ramp(m,j),gy_ramp(m,j));
                    seq_sk.addBlock(gx_spoil,gy_spoil,gz_spoil);
                elseif params.gen.ro_type == 'c'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq.addBlock(gx_pre,gy_pre)
                    else
                        seq.addBlock(gx_pre,gy_pre,tmp_blip)
                    end
                    for k = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
                        gx.amplitude = -gx.amplitude;
                        if k == 1
                            seq.addBlock(gx,gy_blip_up,adc);
                        else
                            seq.addBlock(gx,gy_blips,adc);
                        end
                    end
                    if gx.amplitude > 0
                        gx.amplitude = -gx.amplitude;
                    end
                    seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                    seq.addBlock(tr_delay);
                end
                seq_sk.addBlock(sk_min_tr_delay); 
            end
        end
    end
    seq_sk.addBlock(sk_post_delay);    % Skope post delay
end

%% Actual scan
seq = mr.Sequence();

%%%%% SS-SI-VASO:
if params.gen.seq == 1
% Calibration scan for EPI
%     if params.gen.ro_type == 'c'
%         seq.addBlock(gx_pre)
%         for k = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
%             gx.amplitude = -gx.amplitude;
%             seq.addBlock(gx,adc);
%         end
%         if gx.amplitude > 0
%             gx.amplitude = -gx.amplitude;
%         end
%         seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
%         seq.addBlock(tr_delay);
%     end
    if params.vaso.foci; seq.addBlock(rf_foci); end              % FOCI
    if params.vaso.f_v_delay > 0; seq.addBlock(f_v_delay); end   % FOCI-VASO del%MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslabay
%     if params.gen.fat_sat; seq.addBlock(rf_fs,gz_fs);   end      % fat-sat
    seq.addBlock(ext_trig);                                      % External trigger
    % VASO readout
    for i=1:params.gen.n(3)
        if i==1; tr0 = seq.duration(); end                           % save seq dur to calc TR
        if params.gen.fat_sat; seq.addBlock(rf_fs,gz_fs);   end      % fat-sat
        for j=1:params.gen.seg
            rf = rf0;
            rf.phaseOffset = rf_phase_offset(i);
            adc.phaseOffset = adc_phase_offset(i);
            seq.addBlock(rf,gz);
            if i==1; te0 = seq.duration(); end                           % save seq dur to calc TE
            seq.addBlock(gzReph);
            if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
            % EPI navigators
            if params.gen.ro_type == 'c'
                for i_nav = 1:3
                    gx.amplitude = -gx.amplitude;
                    seq.addBlock(gx,adc); 
                end
            end
            % Gz blip
            if i < floor((params.gen.n(3)/2)+1)
%                     seq.addBlock(gz_blips(i));
                    tmp_blip = gz_blips(i);
            elseif i > floor((params.gen.n(3)/2)+1)
%                     seq.addBlock(gz_blips(i-1)); 
                    tmp_blip = gz_blips(i-1);
            end
            if i==1; te1 = seq.duration(); end                           % save seq dur to calc TE
            if params.gen.ro_type == 's'
                if i == floor((params.gen.n(3)/2)+1)
                    seq.addBlock(gx(i,j),gy(i,j),adc);
                else
                    seq.addBlock(gx(i,j),gy(i,j),tmp_blip,adc);
                end
                seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
            elseif params.gen.ro_type == 'c'
                if i == floor((params.gen.n(3)/2)+1)
                    seq.addBlock(gx_pre,gy_pre)
                else
                    seq.addBlock(gx_pre,gy_pre,tmp_blip)
                end
                for k = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
                    gx.amplitude = -gx.amplitude;
                    if k == 1
                        seq.addBlock(gx,gy_blip_up,adc);
                    else
                        seq.addBlock(gx,gy_blips,adc);
                    end
                end
                if gx.amplitude > 0
                    gx.amplitude = -gx.amplitude;
                end
                seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                seq.addBlock(tr_delay);
            end
    %         seq.addBlock(dummy_delay);          % AMM: Temp: Adding a delay to let mag recover
        end
    if i==1; tr1 = seq.duration(); end                           % save seq dur to calc TR
    end
    if params.vaso.v_b_delay > 0; seq.addBlock(v_b_delay); end   % VASO-BOLD delay
    seq.addBlock(ext_trig);                                      % External trigger
    % BOLD readout
    for i=1:params.gen.n(3)
        if params.gen.fat_sat; seq.addBlock(rf_fs,gz_fs);   end      % fat-sat
        for j=1:params.gen.seg
            rf = rf0;
            rf.phaseOffset = rf_phase_offset(i);
            adc.phaseOffset = adc_phase_offset(i);
            seq.addBlock(rf,gz);
            seq.addBlock(gzReph);
            if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
            % EPI navigators
            if params.gen.ro_type == 'c'
                for i_nav = 1:3
                    gx.amplitude = -gx.amplitude;
                    seq.addBlock(gx,adc); 
                end
            end
            % Gz blip
            if i < floor((params.gen.n(3)/2)+1)
%                     seq.addBlock(gz_blips(i));
                    tmp_blip = gz_blips(i);
            elseif i > floor((params.gen.n(3)/2)+1)
%                     seq.addBlock(gz_blips(i-1)); 
                    tmp_blip = gz_blips(i-1);
            end
            if params.gen.ro_type == 's'
                if i == floor((params.gen.n(3)/2)+1)
                    seq.addBlock(gx(i,j),gy(i,j),adc);
                else
                    seq.addBlock(gx(i,j),gy(i,j),tmp_blip,adc);
                end
                seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
            elseif params.gen.ro_type == 'c'
                if i == floor((params.gen.n(3)/2)+1)
                    seq.addBlock(gx_pre,gy_pre)
                else
                    seq.addBlock(gx_pre,gy_pre,tmp_blip)
                end
                for k = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
                    gx.amplitude = -gx.amplitude;
                    if k == 1
                        seq.addBlock(gx,gy_blip_up,adc);
                    else
                        seq.addBlock(gx,gy_blips,adc);
                    end
                end
                if gx.amplitude > 0
                    gx.amplitude = -gx.amplitude;
                end
                seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                seq.addBlock(tr_delay);
            end
    %         seq.addBlock(dummy_delay);          % AMM: Temp: Adding a delay to let mag recover
        end
    end
    if params.vaso.b_f_delay > 0; seq.addBlock(b_f_delay); end   % BOLD-FOCI delay
%%%%%%  ABC
elseif params.gen.seq == 2
%     % Calibration scan for EPI
%     if params.gen.ro_type == 'c'
%         seq.addBlock(gx_pre)
%         for k = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
%             gx.amplitude = -gx.amplitude;
%             seq.addBlock(gx,adc);
%         end
%         if gx.amplitude > 0
%             gx.amplitude = -gx.amplitude;
%         end
%         seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
%         seq.addBlock(tr_delay);
%     end
    % MT pulse
    if params.mt.mt == 1
        seq.addBlock(MT);                                           % MT pulse 
        seq.addBlock(gz_spoil);
    end
    % if params.gen.fat_sat; seq.addBlock(rf_fs,gz_fs);   end      % fat-sat
    seq.addBlock(ext_trig);                                      % External trigger
    if i==1; tr0 = seq.duration(); end                           % save seq dur to calc TR
    % ABC readout
    dur0 = 0;
    for i=1:params.gen.n(3)
        if params.gen.fat_sat; seq.addBlock(rf_fs,gz_fs);   end      % fat-sat
        for j=1:params.gen.seg
            rf = rf0;
            rf.phaseOffset = rf_phase_offset(i);
            adc.phaseOffset = adc_phase_offset(i);
            seq.addBlock(rf,gz);
            if i==1; te0 = seq.duration(); end                           % save seq dur to calc TE
            seq.addBlock(gzReph);
            if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
            % EPI navigators
            if params.gen.ro_type == 'c'
                for i_nav = 1:3
                    gx.amplitude = -gx.amplitude;
                    seq.addBlock(gx,adc); 
                end
            end
            % Gz blip
            if i < floor((params.gen.n(3)/2)+1)
%                     seq.addBlock(gz_blips(i));
                    tmp_blip = gz_blips(i);
            elseif i > floor((params.gen.n(3)/2)+1)
%                     seq.addBlock(gz_blips(i-1)); 
                    tmp_blip = gz_blips(i-1);
            end
            if i==1; te1 = seq.duration(); end                           % save seq dur to calc TE
            % I am changing this a lot, should be fixed......
            if params.gen.ro_type == 's'
                if params.spi.type == 1
                    % Here I need prewiders.. to take me away from the
                    % center of kspace
                    if i == floor((params.gen.n(3)/2)+1)
                        seq.addBlock(gx_pre,gy_pre)
                        seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                        seq.addBlock(gx(i,j),gy(i,j),adc);
                    else
                        seq.addBlock(gx_pre,gy_pre,tmp_blip)
                        seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                        seq.addBlock(gx(i,j),gy(i,j),adc);
                    end
                elseif params.spi.type == 0
                    if i == floor((params.gen.n(3)/2)+1)
                        seq.addBlock(gx(i,j),gy(i,j),adc);
                    else
                        seq.addBlock(gx(i,j),gy(i,j),tmp_blip,adc);
                    end
                        seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                end
                seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
            elseif params.gen.ro_type == 'c'
                if i == floor((params.gen.n(3)/2)+1)
                    seq.addBlock(gx_pre,gy_pre)
                else
                    seq.addBlock(gx_pre,gy_pre,tmp_blip)
                end
                for k = 1:params.gen.n(2)/params.epi.ry
                    for m = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
                        gx.amplitude = -gx.amplitude;
                        if m == 1
                            seq.addBlock(gx,gy_blip_up,adc);
                        else
                            seq.addBlock(gx,gy_blips,adc);
                        end
                    end
                end
                seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                seq.addBlock(tr_delay);
            end
            if i==1; tr1 = seq.duration(); end                           % save seq dur to calc TR
            % Adding MT pulse every ~180ms
            dur1 = seq.duration();
            if params.mt.mt == 1
                if dur1-dur0 > params.mt.mt_rep && j == params.gen.seg
                    seq.addBlock(MT);                                           % MT pulse 
                    seq.addBlock(gz_spoil);
                    dur0 = dur1;
                end
            end
            % seq.addBlock(dummy_delay);          % AMM: Temp: Adding a delay to let mag recover
        end
    end
%%%%%%  Fieldmap
elseif params.gen.seq == 3
    % Calibration scan for EPI
    if params.gen.ro_type == 'c'
        seq.addBlock(gx_pre)
        for k = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
            gx.amplitude = -gx.amplitude;
            seq.addBlock(gx,adc);
        end
        if gx.amplitude > 0
            gx.amplitude = -gx.amplitude;
        end
        seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
        seq.addBlock(tr_delay);
    end
    % TE loop
    for l = 1:length(params.epi.te)
        % Readout
        for i=1:params.gen.n(3)
            for j=1:params.gen.seg
                rf = rf0;
                rf.phaseOffset = rf_phase_offset(i);
                adc.phaseOffset = adc_phase_offset(i);
                seq.addBlock(rf,gz);
                seq.addBlock(gzReph);
                if params.epi.te > 0; seq.addBlock(fm_te_delay(l)); end
                if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
                % EPI navigators
                if params.gen.ro_type == 'c'
                    for i_nav = 1:3
                        gx.amplitude = -gx.amplitude;
                        seq.addBlock(gx,adc); 
                    end
                end
                % Gz blip
                if i < floor((params.gen.n(3)/2)+1)
    %                     seq.addBlock(gz_blips(i));
                        tmp_blip = gz_blips(i);
                elseif i > floor((params.gen.n(3)/2)+1)
    %                     seq.addBlock(gz_blips(i-1)); 
                        tmp_blip = gz_blips(i-1);
                end
                if params.gen.ro_type == 's'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq.addBlock(gx(i,j),gy(i,j),adc);
                    else
                        seq.addBlock(gx(i,j),gy(i,j),tmp_blip,adc);
                    end
                    seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                    seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                elseif params.gen.ro_type == 'c'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq.addBlock(gx_pre,gy_pre)
                    else
                        seq.addBlock(gx_pre,gy_pre,tmp_blip)
                    end
                    for k = 1:round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
                        gx.amplitude = -gx.amplitude;
                        if k == 1
                            seq.addBlock(gx,gy_blip_up,adc);
                        else
                            seq.addBlock(gx,gy_blips,adc);
                        end
                    end
                    if gx.amplitude > 0
                        gx.amplitude = -gx.amplitude;
                    end
                        seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
                        seq.addBlock(tr_delay);
                end
        %         seq.addBlock(dummy_delay);          % AMM: Temp: Adding a delay to let mag recover
            end
        end
        if params.vaso.b_f_delay > 0; seq.addBlock(b_f_delay); end   % BOLD-FOCI delay  
    end
end

%% check whether the timing of the sequence is correct
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
    plane_samples = adc.numSamples*params.spi.interl;
elseif params.gen.ro_type == 'c'
    plane_samples = adc.numSamples*round(params.gen.n(2)/params.epi.ry*params.epi.pf);
    % Discarding the EPI navigator samples, here I have 3
    j = j+(adc.numSamples*3);
    tmp = ktraj_adc(1,j:j+adc.numSamples*3);
    tmp_mx = max(tmp(:));
end
for i=1:params.gen.n(3)
        ks_traj.kx(:,i) = ktraj_adc(1,j:j+plane_samples-1);
        ks_traj.ky(:,i) = ktraj_adc(2,j:j+plane_samples-1);
        ks_traj.kz(:,i) = ktraj_adc(3,j:j+plane_samples-1);
        j = j+plane_samples;
        if params.gen.ro_type == 'c'
            j = j+(adc.numSamples*3);
        end
end

% Plotting traj
figure(15);
hold on
plot3(ks_traj.kx(:),ks_traj.ky(:),ks_traj.kz(:)); title('3D K-space'); xlabel('Kx (1/m)'); ylabel('Ky (1/m)'); zlabel('Kz (1/m)'); 
view(2)

%% Adding some extra parameters to params
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
tmp = mod(params.gen.n,2); params.gen.n = params.gen.n+tmp;

%% Adding ro_samples, acqTR,volTR, pairTR, TE and other params
params.gen.ro_samples = adc.numSamples;
params.gen.TE = te1-te0+(mr.calcDuration(rf)/2);
params.gen.acqTR = tr1-tr0;
params.gen.volTR = (params.gen.acqTR.*params.gen.n(3));
params.gen.pairTR = params.vaso.f_v_delay+rf_foci.shape_dur+(params.gen.volTR*2);
% params.gen.pairTR = params.gen.volTR*2;
params.gen.adc_dwell = adc.dwell;

%% Calculating time vector for MRIReco.jl
% AMM: Here I might need to take into account the gradient delays...
% julia_time = params.gen.TE+adc.delay+(adc.dwell:adc.dwell:adc.duration);
julia_time = repmat(params.gen.TE+(0:adc.dwell:adc.duration-adc.dwell),1,params.spi.interl);
params.gen.t_vector = julia_time;

%% Designed segment length
seq.setDefinition('MaxAdcSegmentLength',params.gen.adc_split);

%% Saving files
% Check if folder exist, if no it creates it
if exist(sprintf('./data/%s',folder_name)) == 0
    tmp = sprintf('./data/%s',folder_name);
    system(sprintf('mkdir %s',tmp));
    system(sprintf('mkdir %s/acq',tmp));
    system(sprintf('mkdir %s/acq/romeo',tmp));
    system(sprintf('mkdir %s/analysis',tmp));
    system(sprintf('mkdir %s/ismrmd',tmp));
    system(sprintf('mkdir %s/ismrmd/2d',tmp));
    system(sprintf('mkdir %s/ismrmd/3d',tmp));
    system(sprintf('mkdir %s/raw',tmp));
    system(sprintf('mkdir %s/raw/twix',tmp));
    system(sprintf('mkdir %s/raw/dicom',tmp));
    system(sprintf('mkdir %s/recon',tmp));
    system(sprintf('mkdir %s/recon/2d',tmp));
    system(sprintf('mkdir %s/recon/3d',tmp));
    system(sprintf('mkdir %s/tmp',tmp));
end

namestr = strcat('./data/',folder_name,'/acq/',seq_name);
seq.write(strcat(namestr,'.seq'));
if params.gen.skope == 1
    seq_sk.write(strcat(namestr,'_sk.seq'));
end
save(strcat(namestr,'_ks_traj_nom.mat'),'ks_traj')
save(strcat(namestr,'_params.mat'),'params')
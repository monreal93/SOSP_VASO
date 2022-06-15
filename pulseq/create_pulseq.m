clear all; clc
% Change directory to path where sosp_vaso git has been cloned
cd /home/amonreal/Documents/PhD/PhD_2022/sosp_vaso/
addpath(genpath("./pulseq/functions"))
% Add path to pulseq installation folder
addpath(genpath("/home/amonreal/Documents/PhD/tools/pulseq/"))

%% ToDo
% - Check exactly how the saturation pulse should be, what spoliers I need?
% - Check why rfDeadtime in VM has to be at least 150
% - How am I gonna take into account the gradient delays (trajectory,adc
% and time vector)
% - Check if I need a delay for the skope trigger

%% Define parameters
% Set system limits (MaxGrad=70, MaxSlew = 200);
lims = mr.opts('MaxGrad',65,'GradUnit','mT/m',...
    'MaxSlew',190,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6,'rfDeadtime', 180e-6,'adcDeadTime', 10e-6);  % To read it in VM I need rfDeadtime = 180e-6

folder_name = 'sv_06142022';            % Day I am scanning
seq_name = 'sample';                      % use sv_n (n for the diff scans at each day)

% General parameters
params.gen.fov = [200 200 24].*1e-3;
params.gen.res = [0.8 0.8 0.8].*1e-3;
params.gen.fa = 17;
params.gen.te = 0e-3;
params.gen.ro_type = 's';           % 's'-Spiral, 'c'-Cartesian
params.gen.kz = 1;                  % Acceleration in Kz
params.gen.pf = 1;                  % Partial fourier
params.gen.fat_sat = 1;             % Fat saturation (1=yes,0=no)
params.gen.skope = 1;               % Add skope sync scan and triggers

% Spiral parameters
params.spi.type = 0;                % spiral type (for now only spiral out)
params.spi.rotate = 'none';         % Spiral ro0.tation (for now only none)
params.spi.increment = 'linear';    % Spiral increment mode (for now only linear)
params.spi.max_grad  = 65;    % 65     % Peak gradient amplitude for spiral   
params.spi.max_sr = 165;      %165      % Max gradient slew rate for spiral
params.spi.interl = 2;              % Spiral interleaves (for now odata1nly 1)
params.spi.vd = 1.6;                % Variability density
params.spi.rxy = 2.8;                 % In-plane undersampling

% EPI parameters
params.epi.ry = 1;

% VASO parameters
params.vaso.foci = 1;               % FOCI inversion?
params.vaso.tr = 4500e-3;           % volume TR
params.vaso.ti1 = 1800e-3;          % VASO TI1, 
params.vaso.ti2 = params.vaso.ti1+(params.vaso.tr/2);
params.vaso.f_v_delay = 900e-3;     % FOCI-VASO delay
params.vaso.v_b_delay = 0e-3;       % VASO-BOLD delay
params.vaso.b_f_delay = 100e-3;       % BOLD-FOCI delay

% Some calculations
params.gen.del_k = (1./params.gen.fov)*params.gen.kz;
params.gen.n = round(params.gen.fov./params.gen.res);
params.gen.n(3) = params.gen.n(3)/params.gen.kz;
idx = mod(params.gen.n,2)==1;
params.gen.n(idx) = params.gen.n(idx)+1;

%% Create blocks (I will get this into functions)
% TR-FOCI 
[B1_foci,phase,rf_complex,Gz_foci] = tr_foci(3.32,10410e-6);  % Here not sure where this values come from, should I hard code them?
rf_foci = mr.makeArbitraryRf(rf_complex,6.4322, 'system', lims);%, 'Delay',2e-4);

% Fat sat:
% AMM: ToDo: how to properly desing the fat-sat pulse
B0 = 6.98; % 1.5 2.89 3.0
sat_ppm = -3.45;
sat_freq = sat_ppm*1e-6*B0*lims.gamma;
rf_fs = mr.makeGaussPulse(90*pi/180,'system',lims,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
gz_fs = mr.makeTrapezoid('z',lims,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Slice Selective RF gradient 
% ToDo: check values of duration, timeBwProduct
[rf0, gz] = mr.makeSincPulse(params.gen.fa*pi/180,'system',lims,'Duration',2.6e-3,...
    'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2);

%% Preparing readout elements
[spiral_grad_shape,adcSamples,adcDwell,rf_phase_offset,adc_phase_offset,gz_blips_waveform] = prepare_spirals_rf_grad_adc(params,lims);
if params.gen.ro_type == 's'
    % Prepare Spirals
    [spiral_grad_shape,adcSamples,adcDwell,rf_phase_offset,adc_phase_offset,gz_blips_waveform] = prepare_spirals_rf_grad_adc(params,lims);
    for i=1:params.spi.interl
        % Readout gradients
        gx(i) = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:,i), lims);
        gy(i) = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:,i), lims);

        % ADC
        adc = mr.makeAdc(adcSamples,'Dwell',adcDwell);%,'Delay',lims.adcDeadTime);

        % Readout ramp down
        % AMM: need to nicely define the time (0.001) of this ramp gradients
        gx_ramp(i) =mr.makeExtendedTrapezoid('x','times',[0 0.001],'amplitudes',[spiral_grad_shape(1,end,i),0]);
        gy_ramp(i) =mr.makeExtendedTrapezoid('y','times',[0 0.001],'amplitudes',[spiral_grad_shape(2,end,i),0]);
    end
elseif params.gen.ro_type == 'c'
    % Prepare EPI
    % AMM: Need to properly define the dwell time, now I have 2e-6
    gx = mr.makeTrapezoid('x',lims,'riseTime',1e-4,'flatArea',params.gen.n(1)*params.gen.del_k(1),'FlatTime',params.gen.n(1)*2e-6); % Dwell time 2e-6
    gx_pre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2);
    gy_pre = mr.makeTrapezoid('y',lims,'Area',-(params.gen.del_k(2)*params.gen.n(2))/2);
    gy_blip = mr.makeTrapezoid('y',lims,'Area',params.gen.del_k(2).*params.epi.ry);
    
    % ADC
    adcSamples = params.gen.n(1);
    adcDwell = 2e-6;   % AMM: Here i use 2e-6 need to properly define it
    adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',gx.riseTime);
    
end

%% Prepare kz Blips:
for i=1:params.gen.n(3)
    area = -(params.gen.del_k(3)*(params.gen.n(3)/2))+(params.gen.del_k(3)*(i-1));
    dur = ceil(2*sqrt(area/lims.maxSlew)/10e-6)*10e-6;
    if area ~= 0
        gz_blips(i) = mr.makeTrapezoid('z',lims,'Area',area);
    end
end

gz_blips(round(params.gen.n(3)/2)+1) = [];    % Removing the empty blip...

% spoilers
% AMM: Todo: Need to confirm the size (amplitude/area) of this spoiler
gz_spoil=mr.makeTrapezoid('z',lims,'Area',params.gen.del_k(1)*params.gen.n(1)*4);

% Delays, triggers
if params.gen.te > 0; te_delay = mr.makeDelay(params.gen.te-(max(rf0.t)-min(rf0.t))/2); end         % TE delay
if params.vaso.f_v_delay > 0; f_v_delay = mr.makeDelay(params.vaso.f_v_delay);          end         % FOCI-VASO delay
if params.vaso.v_b_delay > 0; v_b_delay = mr.makeDelay(params.vaso.v_b_delay);          end         % VASO-BOLD delay
if params.vaso.b_f_delay > 0; b_f_delay = mr.makeDelay(params.vaso.b_f_delay);          end         % BOLD-FOCI delay
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
ext_trig = mr.makeDigitalOutputPulse('ext1','duration',10e-6,'system',lims);                        % External trigger
dummy_delay = mr.makeDelay(10e-3);           % AMM: Temp: Delay to use as dummy anywhere

%% Add blocks to Skope seq
if params.gen.skope == 1
    seq_sk=mr.Sequence();          % Create a new sequence object
    seq_sk.addBlock(sk_pre_delay); % Skope initial delay
    % Partitions loop
    for i=1:params.gen.n(3)
        % Spiral Interleaves loop
        for j = 1:params.spi.interl
            adc.phaseOffset = adc_phase_offset(i);
            seq_sk.addBlock(gzReph);
            if params.gen.te > 0; seq_sk.addBlock(te_delay); end % TE delay
            seq_sk.addBlock(skope_trig);
            seq_sk.addBlock(sk_int_delay);     % Gradient free interval
            if i < (params.gen.n(3)/2)+1
                seq_sk.addBlock(gz_blips(i));
            elseif i > (params.gen.n(3)/2)+1
                seq_sk.addBlock(gz_blips(i-1)); 
            end
            if params.gen.ro_type == 's'
                seq_sk.addBlock(gx(j),gy(j),adc);
                seq_sk.addBlock(gx_ramp(j),gy_ramp(j),gz_spoil);
            elseif params.gen.ro_type == 'c'
                seq_sk.addBlock(gx_pre,gy_pre)
                for k = 1:params.gen.n(2)/params.epi.ry
                    seq_sk.addBlock(gx,adc)
                    gx.amplitude = -gx.amplitude;
                    seq_sk.addBlock(gy_blip)
                end
                seq_sk.addBlock(gz_spoil);
            end
            seq_sk.addBlock(sk_min_tr_delay); 
        end
    end
    seq_sk.addBlock(sk_post_delay);    % Skope post delay
end

%% Actual scan
seq = mr.Sequence();
if params.vaso.foci; seq.addBlock(rf_foci); end                                          % FOCI
if params.vaso.f_v_delay > 0; seq.addBlock(f_v_delay); end   % FOCI-VASO del%MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslabay
if params.gen.fat_sat; seq.addBlock(rf_fs,gz_fs);   end      % fat-sat  
    % VASO
    for i=1:params.gen.n(3)
        for j=1:params.spi.interl
            rf = rf0;
            rf.phaseOffset = rf_phase_offset(i);
            adc.phaseOffset = adc_phase_offset(i);
            seq.addBlock(rf,gz);
            seq.addBlock(gzReph);
            if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
            if i < (params.gen.n(3)/2)+1
                    seq.addBlock(gz_blips(i));
            elseif i > (params.gen.n(3)/2)+1
                    seq.addBlock(gz_blips(i-1)); 
            end
            if params.gen.ro_type == 's'
                seq.addBlock(gx(j),gy(j),adc);
                seq.addBlock(gx_ramp(j),gy_ramp(j),gz_spoil);
            elseif params.gen.ro_type == 'c'
                seq.addBlock(gx_pre,gy_pre)
                for k = 1:params.gen.n(2)/params.epi.ry
                    seq.addBlock(gx,adc)
                    gx.amplitude = -gx.amplitude;
                    seq.addBlock(gy_blip)
                end
                seq.addBlock(gz_spoil);
            end
    %         seq.addBlock(dummy_delay);          % AMM: Temp: Adding a delay to let mag recover
        end
    end
    if params.vaso.v_b_delay > 0; seq.addBlock(v_b_delay); end   % VASO-BOLD delay
    % BOLD
    for i=1:params.gen.n(3)
        for j=1:params.spi.interl
            rf = rf0;
            rf.phaseOffset = rf_phase_offset(i);
            adc.phaseOffset = adc_phase_offset(i);
            seq.addBlock(rf,gz);
            seq.addBlock(gzReph);
            if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
            if i < (params.gen.n(3)/2)+1
                seq.addBlock(gz_blips(i));
            elseif i > (params.gen.n(3)/2)+1
                seq.addBlock(gz_blips(i-1)); 
            end
            if params.gen.ro_type == 's'
                seq.addBlock(gx(j),gy(j),adc);
                seq.addBlock(gx_ramp(j),gy_ramp(j),gz_spoil);
            elseif params.gen.ro_type == 'c'
                seq.addBlock(gx_pre,gy_pre)
                for k = 1:params.gen.n(2)/params.epi.ry
                    seq.addBlock(gx,adc)
                    gx.amplitude = -gx.amplitude;
                    seq.addBlock(gy_blip)
                end
                seq.addBlock(gz_spoil);
            end
    %         seq.addBlock(dummy_delay);          % AMM: Temp: Adding a delay to let mag recover
        end
    end
    if params.vaso.b_f_delay > 0; seq.addBlock(b_f_delay); end   % BOLD-FOCI delay

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
% [ktraj_adc, ktraj, t_excitation, t_refocusing] = seq.calculateKspace();
% if params.gen.skope
%     j=adc.numSamples*params.gen.n(3)+1;
% else
    j = 1;
% end

if params.gen.ro_type == 's'
    plane_samples = adc.numSamples*params.spi.interl;
elseif params.gen.ro_type == 'c'
    plane_samples = adc.numSamples*params.gen.n(2)/params.epi.ry;
end
for i=1:params.gen.n(3)
        ks_traj.kx(:,i) = ktraj_adc(1,j:j+plane_samples-1);
        ks_traj.ky(:,i) = ktraj_adc(2,j:j+plane_samples-1);
        ks_traj.kz(:,i) = ktraj_adc(3,j:j+plane_samples-1);
        j = j+plane_samples;
end

% Plotting traj
figure(17);
hold on
plot3(ks_traj.kx(:),ks_traj.ky(:),ks_traj.kz(:)); title('3D K-space'); xlabel('Kx (1/m)'); ylabel('Ky (1/m)'); zlabel('Kz (1/m)'); 
view(2)

% seq.plot()

%% Adding some extra parameters to params
params.spi.ro_samples = adc.numSamples;
if params.gen.ro_type == 's'
    params.gen.acqTR = rf.shape_dur+mr.calcDuration(gzReph) ...
        +mr.calcDuration(gz_blips(1))+adc.duration+gx_ramp(1).shape_dur;
    params.gen.TE = (gz.flatTime/2)+gz.fallTime+mr.calcDuration(gzReph)+mr.calcDuration(gz_blips(1));
elseif params.gen.ro_type == 'c'
    params.gen.acqTR = rf.shape_dur+mr.calcDuration(gzReph) ...
        + mr.calcDuration(gz_blips(1))...
        + (mr.calcDuration(gx)+mr.calcDuration(gy_blip))*(params.gen.n(2)/params.epi.ry);
    params.gen.TE = (gz.flatTime/2)+gz.fallTime+mr.calcDuration(gzReph)+mr.calcDuration(gz_blips(1));
end
params.gen.volTR = (params.gen.acqTR.*params.gen.n(3));
params.gen.pairTR = params.vaso.f_v_delay+rf_foci.shape_dur+(params.gen.volTR*2);


%% Calculating time vector for MRIReco.jl
% AMM: Here I might need to take into account the gradient delays...
% julia_time = params.gen.TE+adc.delay+(adc.dwell:adc.dwell:adc.duration);
julia_time = repmat(params.gen.TE+(0:adc.dwell:adc.duration-adc.dwell),1,params.spi.interl);
params.gen.t_vector = julia_time;

%% Saving files
% Check if folder exist, if no it creates it
if exist(sprintf('./data/%s',folder_name)) == 0
    tmp = sprintf('./data/%s',folder_name);
    system(sprintf('mkdir %s',tmp))
    system(sprintf('mkdir %s/acq',tmp))
    system(sprintf('mkdir %s/ismrmd',tmp))
    system(sprintf('mkdir %s/raw',tmp))
    system(sprintf('mkdir %s/recon',tmp))
    system(sprintf('mkdir %s/tmp',tmp))
end

namestr = strcat('./data/',folder_name,'/acq/',seq_name);
seq.write(strcat(namestr,'.seq'));
if params.gen.skope == 1
    seq_sk.write(strcat(namestr,'_sk.seq'));
end
save(strcat(namestr,'_ks_traj_nom.mat'),'ks_traj')
save(strcat(namestr,'_params.mat'),'params')


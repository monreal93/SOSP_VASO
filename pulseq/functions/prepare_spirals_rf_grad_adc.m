function  [spiral_grad_shape,adcSamples,adcDwell,rf_phase_offset,adc_phase_offset] = prepare_spirals_rf_grad_adc(params,lims)
        
        
        
        %% 3D RF and ADC phase
        % This part is coming from ./Emily/RF_design/calcRFphase_3D.m
        % I don't know how this works for 3D; I kinda know now
        delta_x = 0;
        Gslice = 9.5571;
        Tpulse = 2560;
        Asym = 0.5;
        InitPhase = 90;
        
        f = floor(42.5756*delta_x*Gslice+0.5);
            RF.freqOffset = f;
            RF.InitPhaseSet = -(f*360/1e6*Tpulse*Asym)+InitPhase;
            RF.InitPhaseNeg = -(f*360/1e6*Tpulse*(1-Asym))-InitPhase;
            
        % AMM: Original phase...
        for i = 1:params.gen.n(3) %MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslab
            rf_phase(1,i) = RF.InitPhaseSet+25*i*i+175*i+300;
            rf_phase(2,i) = RF.InitPhaseNeg+25*i*i+175*i+300;
            rf_phase_offset(i) = mod(rf_phase(1,i),360)*pi/180;
            
            adc_phase_offset(i) = rf_phase(1,i)-RF.InitPhaseSet;
            adc_phase_offset(i) = mod(adc_phase_offset(i),360)*pi/180;     
        end
        
%         % AMM: Temp: Trying to start the first pahse with i=0 and not i=1
%         for i = 1:params.gen.n(3) %MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslab
%             rf_phase(1,i) = RF.InitPhaseSet+25*(i-1)*(i-1)+175*(i-1)+300;
%             rf_phase(2,i) = RF.InitPhaseNeg+25*(i-1)*(i-1)+175*(i-1)+300;
%             rf_phase_offset(i) = mod(rf_phase(1,i),360)*pi/180;
%             
%             adc_phase_offset(i) = rf_phase(1,i)-RF.InitPhaseSet;
%             adc_phase_offset(i) = mod(adc_phase_offset(i),360)*pi/180;     
%         end
        
        %% Readout Gradients
        % define k-space parameters
    %     deltak =1/params.gen.fov(1);
        kRadius = round(params.gen.n(1)/2);
        kSamples=round(2*pi*kRadius)*2;
%         kSamples = floor(kSamples./1024)*1024;
        % readoutTime = 4.2e-4;

        % AMM Implementation raw Archimiedean spiral
        k_fov = 1./params.gen.res;
%         last_tetha = (2*pi*k_fov(1)*params.gen.fov(1)/2)+(params.spi.rxy*params.spi.vd*2*pi);
        last_tetha = (2*pi*k_fov(1)*params.gen.fov(1)/2);
        
        if params.spi.rxy == 3
%             last_tetha = last_tetha+20;
        end
        
        last_tetha = floor(last_tetha/(4*pi))*4*pi;
%         last_tetha = floor(last_tetha/2*pi)*2*pi;
        
        tetha = 0:2*pi/1e2:last_tetha;
        
        fov_vd = linspace(params.gen.fov(1)*params.spi.vd,params.gen.fov(1),length(tetha));
        ka = (tetha./(2*pi*fov_vd)).*exp(1i.*(tetha./params.spi.rxy./params.spi.interl));
%         ka = ka.';
%         ka = (tetha./(2*pi*fov_vd))*params.spi.rxy.*exp(1i.*tetha);
        % figure; hold on; plot(ka)

        % Partitions loop
        for j=1:params.gen.n(3)
            % Spiral segments
            for i=1:params.spi.interl
                % Getting k-space trajectory
                tmp = (2*pi/params.spi.interl)*i;
                kaa = ka.*exp(1i*tmp);  % Interleave rotations
                
                if contains(params.spi.rotate,'golden')
                    tmp1 = (10*pi/13)*(j-1);
                    kaa = kaa.*exp(1i*tmp1); % Partition rotations
                elseif contains(params.spi.rotate,'linear')
                    tmp1 = 0;
                    kaa = kaa.*exp(1i*tmp1); % Partition rotations
                end

                kaa=[real(kaa); imag(kaa)];
    
    %             %%%%% Calculating Gradients with Lustig approach %%%%%%%
    %             rv = 16; T = 4e-3; ds = -1;
    %             g_max = params.spi.max_grad*10/100;  % convert mT/m -> G/cm
    %             sr_max = params.spi.max_sr*10/100; % convert mT/m/ms -> G/cm/ms   
    %             C = [squeeze(kaa(1,:)).', squeeze(kaa(2,:)).']./100; % times 100 to make it 1/cm
    %             [C,time,g,s,k] = minTimeGradient(C,rv, 0, 0, g_max, sr_max,T,ds,0);
    %             % Interpolating to the initial size of kaa
    % %             tmp = 0:length(g)-1
    %             tmp = linspace(0,time,length(g));
    %             tmp1 = (time*1e-3)./2e-6.*lims.adcRasterTime;
    %             for j=1:3
    %                 g_new(j,:) = interp1(tmp,g(:,j).',linspace(0,time,length(kaa)));
    %             end
    % 
    %             % Trying to get grad to correct units
    %             g_new = g_new*42.58e6/100;    % convert from G/cm -> Hz/m
    %             g_new = rmmissing(g_new,2);
    %             spiral_grad_shape(:,:,i) = g_new;
    %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
                %%%% Calculating Gradients with Pulseq approach %%%%%%%
                [ga, sa]=mr.traj2grad(kaa);
                
                % AMM: Trying to ramp-up for the gradient, now I am just
                % padding some zeros
                tmp = 190;
                rmp_up = [linspace(0,ga(1,1),tmp); linspace(0,ga(2,1),tmp)];
                ga = [rmp_up ga];
                sa = [rmp_up sa];
                kaa = [zeros(2,tmp) kaa];
                
                
    %             Limit analysis
    %             Using the max limits and a safety_margin
                safety_magrin=0.9; % we need that  otherwise we just about violate the slew rate due to the rounding errors
                dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(lims.maxGrad*safety_magrin)*lims.gradRasterTime;
                dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(lims.maxSlew*safety_magrin))*lims.gradRasterTime;
                dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(params.spi.max_grad*lims.gamma*safety_magrin)*lims.gradRasterTime;
                dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(params.spi.max_sr*lims.gamma*safety_magrin))*lims.gradRasterTime;
                dt_smooth=max([dt_gabs;dt_sabs]);
                dt_min=4*lims.gradRasterTime/kSamples; % we want at least 4 points per revolution
                dt_smooth(dt_smooth<dt_min)=dt_min;
                t_smooth=[0 cumsum(dt_smooth,2)];
                kopt_smooth=interp1(t_smooth, kaa', (0:floor(t_smooth(end)/lims.gradRasterTime))*lims.gradRasterTime)';
                [gos, sos]=mr.traj2grad(kopt_smooth);
                kopt_traj(:,:,i) = kopt_smooth;
    
                spiral_grad_shape(:,:,i,j) = gos;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end    
        end

        % AMM: I am padding some zeros so the gradient shape starts from 0
        % Only needed for Pulseq approach...
        tmp = ceil((size(spiral_grad_shape,2)/2)./100)*2*100;
        if tmp == size(spiral_grad_shape,2)
            tmp = size(spiral_grad_shape,2)+100;
        end
        spiral_grad_shape = padarray(spiral_grad_shape,[0 tmp-size(spiral_grad_shape,2) 0],'pre'); 
        % AMM: Here I am just padding some zeros, so grad starts from 0
%         spiral_grad_shape = padarray(spiral_grad_shape,[0 10 0],'pre'); 

        % figure; plot(kopt_smooth(1,:),kopt_smooth(2,:))
        % figure; plot(spiral_grad_shape(1,:))

        %% ADC
        adcDwell=2e-6; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pulseq approach
        % calculate ADC
        % round-down dwell time to 10 ns
        adcTime = lims.gradRasterTime*size(spiral_grad_shape,2);
        % actually it is trickier than that: the (Siemens) interpreter sequence 
        % per default will try to split the trajectory into segments <=1000 samples
        % and every of these segments will have to have duration aligned to the
        % gradient raster time
        adcSamplesPerSegment=2000; % you may need to play with this number to fill the entire trajectory
%         adcSamplesDesired=kRadius*kSamples; %AMM: what is my actual desired samples?
        adcSamplesDesired=round(adcTime./adcDwell);            % AMM: I am manually setting this params
        adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
        adcSamples=adcSegments*adcSamplesPerSegment;
%         adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
%         adcDwell=2e-6;                          % AMM: Manually setting adcDwell to have 5samples per GradientRaster
        adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 100 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
        if mod(adcSegmentDuration, lims.gradRasterTime)>eps 
            error('ADC segmentation model results in incorrect segment duration');
        end
        % update segment count
%         adcSegments=floor(adcTime/adcSegmentDuration);
        adcSamples=adcSegments*adcSamplesPerSegment;
        
%         % AMM: Another way to get adc Samples
% %         adcSamples = ceil(round(adcTime./adcDwell)/10)*10;
        adcSamples = round(adcTime./adcDwell);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%% Lusing approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         adcSamples = round(time*1e-3./adcDwell);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % AMM: Compensating for the 6 zeros I padded before
%         adcSamples = adcSamples-6;
%         adcSamples = floor(adcSamples/4)*4; % on Siemens the number of ADC samples need to be divisible by 4

        % extend spiral_grad_shape by repeating the last sample
        % this is needed to accomodate for the ADC tuning delay
        spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end,:,:)];
        
        % SOSP Blips
        nz_eff = params.gen.n(3);       % AMM: this should be changed for undersampling kz
        delta_kz = 1/params.gen.fov(3)/1000*params.gen.n(3)/nz_eff;
        kzfov_half = 1/(params.gen.fov(3)/params.gen.n(3))/2/1000; %rad/2pi/mm
        kz_coordinates = (-1*kzfov_half)+delta_kz:delta_kz:kzfov_half;
%         gz_blips = zeros(1,54);
%         for i=1:params.gen.n(3)
%             [gz_blips(i).blip(:),gz_blips(i).amp(:)] = DesignSOSPBlips_AMM(abs(1/kz_coordinates(i)),params);
%             gz_blips1(i).blip(:) = gz_blips(i).blip(:).*gz_blips(i).amp;
%         end


end
function  [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc(params,lims)
        
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
        
        last_tetha = ceil(last_tetha/(2*pi))*2*pi;
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
    
                %%%%% Calculating Gradients with Lustig approach %%%%%%%
                rv = 16; T = 4e-3; ds = -1;
                g_max = params.spi.max_grad*10/100;  % convert mT/m -> G/cm
                sr_max = params.spi.max_sr*10/100; % convert mT/m/ms -> G/cm/ms   
                C = [squeeze(kaa(1,:)).', squeeze(kaa(2,:)).']./100; % times 100 to make it 1/cm

                if params.spi.type == 1
                    C = flip(C,1);
                end

                if j==1 || contains('none',params.spi.rotate) == 0
                    %%%%%%% Unsafe spirals
                    [C,time,g,s,k] = minTimeGradient(C,rv, C(1,1), 0, g_max, sr_max,T,ds,0);
                    % Interpolating to the initial size of kaa
        %             tmp = 0:length(g)-1
                    time = round(time*1e-3,3);
                    tmp = linspace(0,time,length(g));
                    tmp1 = (time*1e-3)./2e-6.*lims.adcRasterTime;
                    for k=1:2
                        % g_new(k,:) = interp1(tmp,g(:,k).',linspace(0,time,length(kaa)));
                        g_new(k,:) = interp1(tmp,g(:,k).',(0:floor(time(end)/lims.gradRasterTime))*lims.gradRasterTime);
                    end
                                    % Trying to get grad to correct units
                    g_new = g_new*42.58e6/100;    % convert from G/cm -> Hz/m
                    g_new = rmmissing(g_new,2);
                    %%%%%%%%%%%%
    
%                     %%%%%% Safe spirals
%                     [g_new(1,:),g_new(2,:),time] = fn_safe_spirals(C,params);
%                     %%%%%%%%%
                end

                spiral_grad_shape(:,:,i,j) = g_new(1:2,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 
%                 %%%% Calculating Gradients with Pulseq approach %%%%%%%
%                 [ga, sa]=mr.traj2grad(kaa);
%                 
%                 % AMM: Trying to ramp-up for the gradient, now I am just
%                 % padding some zeros
%                 tmp = 190;
%                 rmp_up = [linspace(0,ga(1,1),tmp); linspace(0,ga(2,1),tmp)];
%                 ga = [rmp_up ga];
%                 sa = [rmp_up sa];
%                 kaa = [zeros(2,tmp) kaa];
%                 
%                 
%     %             Limit analysis
%     %             Using the max limits and a safety_margin
%                 safety_magrin=0.9; % we need that  otherwise we just about violate the slew rate due to the rounding errors
%                 dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(lims.maxGrad*safety_magrin)*lims.gradRasterTime;
%                 dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(lims.maxSlew*safety_magrin))*lims.gradRasterTime;
%                 dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(params.spi.max_grad*lims.gamma*safety_magrin)*lims.gradRasterTime;
%                 dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(params.spi.max_sr*lims.gamma*safety_magrin))*lims.gradRasterTime;
%                 dt_smooth=max([dt_gabs;dt_sabs]);
%                 dt_min=4*lims.gradRasterTime/kSamples; % we want at least 4 points per revolution
%                 dt_smooth(dt_smooth<dt_min)=dt_min;
%                 t_smooth=[0 cumsum(dt_smooth,2)];
%                 kopt_smooth=interp1(t_smooth, kaa', (0:floor(t_smooth(end)/lims.gradRasterTime))*lims.gradRasterTime)';
%                 [gos, sos]=mr.traj2grad(kopt_smooth);
%                 kopt_traj(:,:,i) = kopt_smooth;
%     
%                 spiral_grad_shape(:,:,i,j) = gos;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end    
        end

        % Padding some zeros so the gradient shape starts from 0 
        if params.spi.type == 0
            tmp = ceil((size(spiral_grad_shape,2)/2)./100)*2*100;
            if tmp == size(spiral_grad_shape,2)
                tmp = size(spiral_grad_shape,2)+100;
            end
            spiral_grad_shape = padarray(spiral_grad_shape,[0 tmp-size(spiral_grad_shape,2) 0],'pre'); 
        end

        if params.spi.type == 1
           spiral_grad_shape = padarray(spiral_grad_shape,[0 10 0],'post');  
        end
        

        % AMM: Here I am just padding some zeros, so grad starts from 0
%         spiral_grad_shape = padarray(spiral_grad_shape,[0 10 0],'pre'); 

        % figure; plot(kopt_smooth(1,:),kopt_smooth(2,:))
        % figure; plot(spiral_grad_shape(1,:))

        %% ADC
        % Round-down dwell time to adcRasterTime
        adcDwell = floor((1./(params.spi.bw)/lims.adcRasterTime))*lims.adcRasterTime;
        % updating Spiral BW
        params.spi.bw = 1./adcDwell;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pulseq approach
%         % calculate ADC
%         % round-down dwell time to 10 ns
%         adcTime = lims.gradRasterTime*size(spiral_grad_shape,2);
%         % actually it is trickier than that: the (Siemens) interpreter sequence 
%         % per default will try to split the trajectory into segments <=1000 samples
%         % and every of these segments will have to have duration aligned to the
%         % gradient raster time
%         adcSamplesPerSegment=2000; % you may need to play with this number to fill the entire trajectory
% %         adcSamplesDesired=kRadius*kSamples; %AMM: what is my actual desired samples?
%         adcSamplesDesired=round(adcTime./adcDwell);            % AMM: I am manually setting this params
%         adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
%         adcSamples=adcSegments*adcSamplesPerSegment;
% %         adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
% %         adcDwell=2e-6;                          % AMM: Manually setting adcDwell to have 5samples per GradientRaster
%         adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 100 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
%         if mod(adcSegmentDuration, lims.gradRasterTime)>eps 
%             error('ADC segmentation model results in incorrect segment duration');
%         end
%         % update segment count
% %         adcSegments=floor(adcTime/adcSegmentDuration);
%         adcSamples=adcSegments*adcSamplesPerSegment;
%         
% % %         % AMM: Another way to get adc Samples
% % % %         adcSamples = ceil(round(adcTime./adcDwell)/10)*10;
% %         adcSamples = round(adcTime./adcDwell);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Lusing approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        adcSamples = round(time./adcDwell);
        % In SIEMENS number of ADC samples should be divisible by 4
        adcSamples = floor(adcSamples/4)*4;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Getting the correct number to split the ADC
        % The ADC obj has to be splitted into N equal parts, with duration multiple
        % of 10us
        for i = 1:50
            if mod(adcSamples,i) == 0 && mod(adcSamples/i*adcDwell,10e-9) == 0 && (adcSamples/i) < 8192 && mod(adcSamples/i,4) == 0
                adcSplit = adcSamples/i;
                break
            end
        end
        params.gen.adc_split = adcSplit;
        params.gen.adc_segments = i;
        
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

    %% Checking forbidden frequencies
    check_forbbiden_fq(squeeze(spiral_grad_shape(1,:,1,1)),false)
    title('Forbidden Frequencies Gx')
    check_forbbiden_fq(squeeze(spiral_grad_shape(2,:,1,1)),false)
    title('Forbidden Frequencies Gy')


end
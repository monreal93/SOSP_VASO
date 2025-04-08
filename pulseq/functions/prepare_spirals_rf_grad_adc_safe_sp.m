function  [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc_safe_sp(params)
        
lims = params.gen.lims;
%% Readout Gradients
% AMM Implementation raw Archimiedean spiral
k_fov = 1./params.gen.res;
last_tetha = (2*pi*k_fov(1)*params.gen.fov(1)/2);

% to match resolution as Cartesian (square k-space) they must measure a
%     diameter of 2/sqrt(pi) ~ 1.13 larger than conventional k-space limits (so
%     that the area of the circle equals the area of the square).
% last_tetha = last_tetha * 1.13;

last_tetha = floor(last_tetha/(2*pi))*2*pi;
tetha = 0:2*pi/1e2:last_tetha;

% Positive of negative VD
if params.spi.vd > 0
    fov_vd = linspace(params.gen.fov(1)*params.spi.vd,params.gen.fov(1),length(tetha));
else
    fov_vd = linspace(params.gen.fov(1),params.gen.fov(1)*params.spi.vd*-1,length(tetha));
end

if params.spi.type == 4
    ka = (tetha./(2*pi*fov_vd)).*exp(1i.*(tetha./(params.spi.rxy*2)./params.spi.interl));
else
    ka = (tetha./(2*pi*fov_vd)).*exp(1i.*(tetha./params.spi.rxy./params.spi.interl));
end

% If IN-OUT with shift between, loop again in interleaves...
if params.spi.type == 4
    in_out_rot = 2;
else
    in_out_rot = 1;
end

% Partitions loop
for j=1:params.gen.n(3)
    % Spiral segments
    for i=1:params.spi.interl*in_out_rot
        % Getting k-space trajectory
        tmp = (2*pi/params.spi.interl)*i;
        % Interleave rotations
        kaa = ka.*exp(1i*tmp); 
        
        % If spiral IN/OUT, rotate the OUT part
        % ToDo: Properly find how when to enter this if... 
        if  params.spi.type == 4 && i == params.spi.interl*in_out_rot
            tmp = 2*pi;
            kaa = kaa.*exp(1i*tmp); 
        end
        
        % Partition rotations
        if contains(params.spi.rotate,'golden')
            tmp1 = (10*pi/13)*(j-1);
            kaa = kaa.*exp(1i*tmp1);
        elseif contains(params.spi.rotate,'linear')
            tmp1 = 0;
            kaa = kaa.*exp(1i*tmp1);
        elseif contains(params.spi.rotate,'90')
            tmp1 = (90*pi/180)*(j-1);
            kaa = kaa.*exp(1i*tmp1);
        elseif contains(params.spi.rotate,'120')
            tmp1 = (120*pi/180)*(j-1);
            kaa = kaa.*exp(1i*tmp1);
        elseif contains(params.spi.rotate,'180')
            tmp1 = pi*(j-1);
            kaa = kaa.*exp(1i*tmp1);
        end

        kaa=[real(kaa); imag(kaa)];

        %%%%% Calculating Gradients with Lustig approach %%%%%%%
%         rv = 16; T = 4e-3; ds = -1; % Original for mex function
        rv = []; T = 4e-3; ds = []; % For Matlab function.
        g_max = params.spi.max_grad*10/100;  % convert mT/m -> G/cm
        sr_max = params.spi.max_sr*10/100; % convert mT/m/ms -> G/cm/ms   
        C = [squeeze(kaa(1,:)).', squeeze(kaa(2,:)).']./100; % times 100 to make it 1/cm

        C(:,3) = zeros(size(C,1),1);  % For Matlab function...

        if j==1 || contains('none',params.spi.rotate) == 0
            % Make Gradient within SR and Grad limits
            [C,time,g,s,k] = minTimeGradient(C,rv, 0, 0, g_max, sr_max,T,ds,0);


            %%%%% Temp: trying to make safe spirals...
            gx = g(:,1);
            gy = g(:,2);
            gx = gx.*1e-2;
            gy = gy.*1e-2;
            kx = cumsum(gx);
            ky = cumsum(gy);
%             kx = C(:,1);
%             ky = C(:,2);
            r = sqrt((kx.^2) + (ky.^2));
            
            ff=(1./(r.*1e-3));
%             ff=(norm(gx)./(r*1e-1));
%             figure(80); hold on; plot(ff); ylim([0,4000]); title('Instantaneous frequency')
%             figure(80); hold on; plot(1:length(kx),repmat(500,[1, length(kx)]),'red')
%             figure(80); hold on; plot(1:length(kx),repmat(600,[1, length(kx)]),'red')
%             figure(80); hold on; plot(1:length(kx),repmat(950,[1, length(kx)]),'red')
%             figure(80); hold on; plot(1:length(kx),repmat(1250,[1, length(kx)]),'red')

            % Finding the points when we enter the forbidden frequencies...
%             fb_fq(1) = length(ff);
            tmp = find(ff>500); fb_fq=tmp(end); % 500
            tmp = find(ff>600); fb_fq(2)=tmp(end); % 600
            tmp = find(ff>950); fb_fq(3)=tmp(end); % 950
            tmp = find(ff>1250); fb_fq(4)=tmp(end); % 1250

%             % Temp.. manually setting them...
%             fb_fq(1) = 3339;
%             fb_fq(2) = 2671;
%             fb_fq(3) = 1378;
%             fb_fq(4) = 908;

            fb_fq = flip(fb_fq);
%             fb_fq = [1,fb_fq];

           %%%%%%%
%             fb_fq(4) = fb_fq(4)+200;

            tmp1 = g(:,1)>-0.05;
            tmp2 = g(:,1)<0.05;
            xx = tmp1.*tmp2;

            yy = zeros(size(xx));

            % rg might need to be adjusted to properly find the values
            rg = 180;
            for i_rg = 1:length(fb_fq)
                    yy(fb_fq(i_rg)-(rg/2):fb_fq(i_rg)+(rg/2),1) = 1;
            end
   
            zz = xx.*yy;
            zz = find(zz);
            zzz = zz(diff(zz)>rg);
            zzz = [zzz;zz(end)];
            zzz = [1;zzz];
            zzz = [zzz;length(ff)];

            fb_fq = zzz.';  % Original
            fb_fq(2:end) = fb_fq(2:end)-1;
%             fb_fq(2) = fb_fq(2) - 88;
% 
%             figure; plot(tmp1.*tmp2)
%             hold on
%             plot(fb_fq,1,'x','Color','red')
%             plot(zzz,1,'o','Color','g')

            g_safe = zeros([1,3]);
            C_orig = C;
            g_orig = g;
            time_orig = time;
            time_safe = 0;
            %%%%%% Gradient scaling....
            % The ones I want to sacale are positions 3 and 5
            sr_scaling = [1,2,2,2,2];
            g_scaling = [1,1.75,1.75,1.75,1.75];
            %%%%%%%%%%%%%
            for i_seg=1:length(fb_fq)-1
                if mod(i_seg,2)
%                     [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg):fb_fq(i_seg+1),:),rv, g_safe(end,1), g_max, g_max, sr_max,T,ds,0);
                    if i_seg == 1
%                         [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg):fb_fq(i_seg+1),:),rv, 0, 0, g_max/scaling(i_seg), sr_max/scaling(i_seg),T,ds,0); % Orig
                        [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg):fb_fq(i_seg+1),:),rv, 0, [],  g_max/g_scaling(i_seg), sr_max/sr_scaling(i_seg),T,ds,0);
                    else
%                         [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg):fb_fq(i_seg+1),:),rv, 0, 0, g_max/scaling(i_seg), sr_max/scaling(i_seg),T,ds,0); % Orig
                          [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg):fb_fq(i_seg+1),:),rv, abs(sqrt((g(end,1)^2)+g(end,2)^2)), [],  g_max/g_scaling(i_seg), sr_max/sr_scaling(i_seg),T,ds,0);
                    end
%                     g = g.*-1;
%                     tmp = find(g(:,1)>g_safe(end,1));
%                     g = g(tmp(1):end,:);
                else
%                     [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg)+1:fb_fq(i_seg+1),:),rv, g_safe(end-1,1), g_max/1.7, g_max/1.7, sr_max,T,ds,0);
%                     [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg):fb_fq(i_seg+1),:),rv, 0, 0, g_max/scaling(i_seg), sr_max/scaling(i_seg),T,ds,0);  % Orig
                      [C,time,g,s,k] = minTimeGradient(C_orig(fb_fq(i_seg):fb_fq(i_seg+1),:),rv,abs(sqrt((g(end,1)^2)+g(end,2)^2)), [],  g_max/g_scaling(i_seg), sr_max/sr_scaling(i_seg),T,ds,0);
%                     g = g.*-1;
%                     tmp = find(g(:,1)>g_safe(end,1));
%                     g = g(tmp(1):end,:);
                end
                

                    %%%%% trying to make a smooth transition...
%                     ggg = [];
%                     for i_dims = 1:2
%                         gg=g(:,i_dims);
%                         samples=24;
%                         tmp = linspace(g(end-(samples-1),i_dims),0,samples/2).';
%                         gg=gg(1:end-samples);
%                         ggg(:,i_dims) = [gg; tmp];
%                     end
%                     ggg(:,3) = zeros(length(ggg),1);
%                 tmp1 = g(:,2)>-0.03;
%                 tmp2 = g(:,2)<0.03;
%                 tmp3 = find(tmp1.*tmp2);
%                 tmp_st = tmp3(1);
% 
%                 tmp1 = g(:,2)>-0.03;
%                 tmp2 = g(:,2)<0.03;
%                 tmp3 = find(tmp1.*tmp2);
%                 tmp_end = tmp3(end)-1;

%                 if i_seg == 5
%                     g_safe = [g_safe; g(tmp_st:end,:)];
%                 else
%                     g_safe = [g_safe; g(tmp_st:tmp_end,:)];
%                 end
%               
%             if i_seg > 1
%                 g_tmp = g(40:end,2);
%                 tmp = find(g_safe(:,2)>g_tmp(1));
%                 tmp = tmp(end);
%                 g_tmp = [flip(g_safe(tmp:end,2));g_tmp(2:end)];
%                 g = g_tmp;
%                 
%                 tmp = linspace(-2.2,0,20);
%                 tmp = tmp.^2;
% 
%             end
            %%%%%%

            g_safe = [g_safe;g];
            time_safe = time_safe+time;

            end

            % Checking the safe spirals..
            
            gx = g_safe(:,1);
            gy = g_safe(:,2);
            gx = gx.*1e-2;
            gy = gy.*1e-2;
            kx = cumsum(gx);
            ky = cumsum(gy);
%             kx = C(:,1);
%             ky = C(:,2);
            r = sqrt((kx.^2) + (ky.^2));
            
            ff=(1./(r.*1e-3));
%             ff=(norm(gx)./(r*1e-1));
%             figure(80); hold on; plot(ff);


            g = g_safe;
            time = time_safe;
%             g = g_orig;
%             time = time_orig;

            %%%%%%

            % Adding some zeros to force gradient to start at 0
            tmp = [g(1,:)/2; g(1,:)/4];
            tmp1 = [g(end,:)/2; g(end,:)/4];
            if params.spi.type == 3
                tmp = flip(tmp);
                tmp = tmp(1,:);
                g = [ g; zeros(1,3)];
            elseif params.spi.type == 4
                if i > 1
                    tmp = g1(:,end)/42.58e6*100;
                    tmp = [tmp; 0]';
                    g(1,:) = tmp;
                end
            else
                g = [zeros(20,3); tmp ; g; tmp1; zeros(20,3)];
            end

            time = round(time*1e-3,3);
            tmp = linspace(0,time,length(g));
            tmp1 = linspace(0,time,ceil(time(end)/lims.gradRasterTime));

            for k=1:2
                g_new(k,:) = interp1(tmp,g(:,k).',tmp1,'spline','extrap');
            end

            % Convert gradient units
            g_new = g_new*42.58e6/100;    % convert from G/cm -> Hz/m
            g_new = rmmissing(g_new,2);

            if params.spi.type == 1
                g_new = flip(g_new,2);
            end

            if params.spi.type == 3
                g1 = [flip(g_new,2).*-1,g_new];
                spiral_grad_shape(:,:,i,j) = g1(1:2,:);
            elseif params.spi.type == 4
                if i == 1
                    g1 = [flip(g_new,2)];
                else
                    g1 = [g1,g_new];
                    spiral_grad_shape(:,:,i-1,j) = g1(1:2,:);
                end
            else
                g1 = g_new;
                spiral_grad_shape(:,:,i,j) = g1(1:2,:);
            end

        elseif params.spi.interl > 1 || contains('none',params.spi.rotate) == 1
            if params.spi.type == 4
                spiral_grad_shape(:,:,1,j) = spiral_grad_shape(:,:,1,1);
            else
                spiral_grad_shape(:,:,i,j) = spiral_grad_shape(:,:,i,1);
            end
        end
    end    
end

% Padding some zeros to make it a nice number
tmp = ceil((size(spiral_grad_shape,2)/2))*2;
tmp = ceil(tmp/100)*100;

if params.spi.type == 0
    spiral_grad_shape = padarray(spiral_grad_shape,[0 tmp-size(spiral_grad_shape,2) 0],'pre'); 
elseif params.spi.type == 1
    spiral_grad_shape = padarray(spiral_grad_shape,[0 tmp-size(spiral_grad_shape,2) 0],'post'); 
end

% Adjusting dwell time...
if params.spi.bw == 0
    % Nyquist
    adcDwell = (1/(params.gen.fov(1)*lims.gamma*params.spi.max_grad*1e-3)); % Nyquist limit
    adcDwell = adcDwell*params.spi.rxy_az;                                  % Nyquist undersampled
else
    % Custom
    adcDwell = 1/params.spi.bw; %custom BW
end

adcDwell = floor(adcDwell/lims.adcRasterTime/2)*lims.adcRasterTime*2;

time_new = size(spiral_grad_shape,2).*lims.adcRasterTime/1e-2;
% gradient should finish at gradRaster time
time_new = round(time_new/lims.gradRasterTime)*lims.gradRasterTime;
adcSamples = round(time_new./adcDwell);
% In SIEMENS number of ADC samples should be divisible by 4
adcSamples = ceil(ceil(adcSamples/4)*4/1000)*1000;

% Checking forbidden frequencies, simple approach
scanner = params.gen.field_strength;
check_forbbiden_fq(squeeze(spiral_grad_shape(1,:,1,1)),scanner,true);
check_forbbiden_fq(squeeze(spiral_grad_shape(2,:,1,1)),scanner,false);

% figure(81); hold on; plot(spiral_grad_shape(1,:,1,1)); title("Gx safe and UNsafe")
% figure(82); hold on; plot(spiral_grad_shape(2,:,1,1)); title("Gy safe and UNsafe")
% figure; plot(cumsum(squeeze(spiral_grad_shape(1,:,1,1))),cumsum(squeeze(spiral_grad_shape(2,:,1,1)))); title("K-space safe spiral")

% 
% %%%%% Temp: trying to make safe spirals...
% gx = spiral_grad_shape(1,:,1,1);
% gy = spiral_grad_shape(2,:,1,1);
% gx = gx.*42.58e-9;
% gy = gy.*42.58e-9;
% kx = cumsum(gx);
% ky = cumsum(gy);
% r = sqrt((kx.^2) + (ky.^2));
% 
% ff=(1./(r.*1e-3));
% ff=ff(10:end);
% figure(80); hold on; plot(ff); ylim([0,4000]);
% figure(80); hold on; plot(1:6000,repmat(500,[1, 6000]),'red')
% figure(80); hold on; plot(1:6000,repmat(600,[1, 6000]),'red')
% figure(80); hold on; plot(1:6000,repmat(950,[1, 6000]),'red')
% figure(80); hold on; plot(1:6000,repmat(1250,[1, 6000]),'red')
% check_forbbiden_fq(squeeze(spiral_grad_shape(1,:,1,1)),scanner,true);


end
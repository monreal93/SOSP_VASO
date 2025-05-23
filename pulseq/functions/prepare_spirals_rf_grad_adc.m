function  [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc(params)
        
lims = params.gen.lims;
%% Readout Gradients
% AMM Implementation raw Archimiedean spiral
k_fov = 1./params.gen.res;
last_tetha_x = (2*pi*k_fov(1)*params.gen.fov(1)/2);
last_tetha_y = (2*pi*k_fov(2)*params.gen.fov(2)/2);

% to match resolution as Cartesian (square k-space) they must measure a
%     diameter of 2/sqrt(pi) ~ 1.13 larger than conventional k-space limits (so
%     that the area of the circle equals the area of the square).
% last_tetha_x = last_tetha_x * 1.13;
% last_tetha_y = last_tetha_y * 1.13;

last_tetha_x = floor(last_tetha_x/(2*pi))*2*pi;
last_tetha_y = floor(last_tetha_y/(2*pi))*2*pi;
tetha_x = 0:2*pi/1e2:last_tetha_x;
tetha_y = 0:2*pi/1e2:last_tetha_y;
% tetha_y = linspace(0,last_tetha_y,length(tetha_x));

% Positive of negative VD
if params.spi.vd > 0
    fov_vd_x = linspace(params.gen.fov(1)*params.spi.vd,params.gen.fov(1),length(tetha_x));
    fov_vd_y = linspace(params.gen.fov(2)*params.spi.vd,params.gen.fov(2),length(tetha_y));
else
    fov_vd_x = linspace(params.gen.fov(1),params.gen.fov(1)*params.spi.vd*-1,length(tetha_x));
    fov_vd_y = linspace(params.gen.fov(2),params.gen.fov(2)*params.spi.vd*-1,length(tetha_y));
end

if params.spi.type == 4
    ka = (tetha_x./(2*pi*fov_vd_x)).*exp(1i.*(tetha_x./(params.spi.rxy*2)./params.spi.interl));
%     ka_x = (tetha_x./(2*pi*fov_vd_x)).*cos((tetha_x./(params.spi.rxy*2)./params.spi.interl));
%     ka_y = (tetha_y./(2*pi*fov_vd_y)).*cos((tetha_y./(params.spi.rxy*2)./params.spi.interl));
else
    ka = (tetha_x./(2*pi*fov_vd_x)).*exp(1i.*(tetha_x./params.spi.rxy./params.spi.interl));
%     ka_x = (tetha_x./(2*pi*fov_vd_x)).*sin((tetha_x./params.spi.rxy./params.spi.interl));
%     ka_y = (tetha_y./(2*pi*fov_vd_y)).*sin((tetha_y./params.spi.rxy./params.spi.interl));
end
% ka = ka_x + (ka_y.*1i);

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
        rv = 16; T = 4e-3; ds = -1; % Original
%         rv = []; T = 4e-3; ds = []; % For Matlab function.
        g_max = params.spi.max_grad*10/100;  % convert mT/m -> G/cm
        sr_max = params.spi.max_sr*10/100; % convert mT/m/ms -> G/cm/ms   
        C = [squeeze(kaa(1,:)).', squeeze(kaa(2,:)).']./100; % times 100 to make it 1/cm

        C(:,3) = zeros(size(C,1),1);  % For Matlab function...

        if j==1 || contains('none',params.spi.rotate) == 0
            % Make Gradient within SR and Grad limits
            [C,time,g,s,k] = minTimeGradient(C,rv, C(1,1), 0, g_max, sr_max,T,ds,0);


%             %%%%%% Temp: trying to make safe spirals
%             tmp = k(:,1) + (k(:,2).*1i);
%             phi = unwrap(angle(tmp));
%             rc = (1/(2.*pi*params.gen.fov(1))).*((((phi.^2)+1).^(3/2))./((phi.^2)+2));
%             rc = rc.';
%             xx = (lims.gamma.*params.spi.max_grad*1e-3./rc(1:end-1));
%             tt = 0:T:time-T;
%             %%%%%%

            % Adding some zeros to force gradient to start at 0
            tmp = [g(1,:)/2; g(1,:)/4];
            tmp1 = [g(end,:)/2; g(end,:)/4];
            if params.spi.type == 3
                tmp = flip(tmp);
                tmp = tmp(1,:);
                g = [ zeros(1,3); g; zeros(1,3)];
            elseif params.spi.type == 4
                if i > 1
                    tmp = g1(:,end)/42.58e6*100;
                    tmp = [tmp; 0]';
                    g(1,:) = tmp;
                else
                    g = padarray(g,[4,0]);
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
time_new = round(time_new/lims.gradRasterTime)*lims.gradRasterTime; % original
time_new = round(time_new/lims.adcRasterTime)*lims.adcRasterTime; % New with adc Raster
adcSamples = round(time_new./adcDwell);
% AMM: ToDo: This makes the ADC unecessary longer, but seems to be needed
% to run at the scanner, fix it, should it be a multiple of 1000? or 100?
% In SIEMENS number of ADC samples should be divisible by 4
% adcSamples = ceil(ceil(adcSamples/4)*4/1000)*1000; % original
% adcSamples = round(round(adcSamples/4)*4/100)*100; % New, rounding to 100
% adcSamples = round(round(adcSamples/4)*4/10)*10; % New, rounding to 10
% adcSamples = round(round(adcSamples/4)*4); % New, No rounding to multiple of 10
if params.spi.type ~= 3 % Don't do it for spiral IN-OUT
    adcSamples = round(round(adcSamples/4)*4/200)*200; % New, rounding to 200
end

% Checking forbidden frequencies, simple approach
scanner = params.gen.field_strength;
check_forbbiden_fq(squeeze(spiral_grad_shape(1,:,1,1)),scanner,true);
check_forbbiden_fq(squeeze(spiral_grad_shape(2,:,1,1)),scanner,false);

end
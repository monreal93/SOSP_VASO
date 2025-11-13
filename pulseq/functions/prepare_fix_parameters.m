function [params,ro_blocks] = prepare_fix_parameters(params)

    % Saving nominal resolution before it gets uptaded
    params.gen.nom_res = params.gen.res;
    
    % to match resolution as Cartesian (square k-space) they must measure a
    %     diameter of 2/sqrt(pi) ~ 1.13 larger than conventional k-space limits (so
    %     that the area of the circle equals the area of the square).
    % Reference: Pipe JG, Zwart NR (2014) Spiral trajectory design: A flexible numerical algorithm and base analytical equations. 
    params.gen.res(1:2) = params.gen.res(1:2)./1.13;

    % If ME-GRE.. set params.gen.echos to # of echos/contrasts
    if params.gen.seq == 3
        params.gen.echos = length(params.gen.multi_te);
    end

    params.gen.del_k = (1./params.gen.fov);
    if params.gen.ro_type == 's'; params.gen.del_k(1:2) = params.gen.del_k(1:2).*params.spi.rxy; end
    if params.gen.ro_type == 'c'; params.gen.me_gre =0; end % ME-GRE only for spiral
    if params.gen.ro_type == 'c'; params.spi.interl = 1; end
    if params.gen.ro_type == 'c'; params.gen.del_k(2) = params.gen.del_k(2).*params.epi.ry; end
    if params.mt.bold; params.mt.mt = 1; end        % If MT BOLD corrected, use MT 
    if params.gen.kz == 1; params.gen.kz_caipi = 1; end % If no Rz, no CAIPI kz
    params.gen.del_k(3) = params.gen.del_k(3).*params.gen.kz;
%     if params.gen.kz_caipi == 1; params.gen.del_k(3) = params.gen.del_k(3)/2; end  % If Kz caipi, kz/blip /2
    params.gen.n = round(params.gen.fov./params.gen.res);
    if and(params.gen.seq == 1,params.vaso.bold_ref == 1) || and(params.gen.seq == 2,params.mt.bold)
        ro_blocks = 2; 
%     elseif params.gen.seq == 3
%         ro_blocks = length(params.gen.multi_te);
    else
        ro_blocks = 1; 
    end
    % Partition or phase oversampling
    params.gen.n_ov = params.gen.n;
    params.gen.fov_ov = params.gen.fov;
    if params.gen.ph_oversampling > 0
        % % Increase FOV? Actually increasing FOV...
        % % params.gen.fov(3) = params.gen.fov(3).*(1+(params.gen.ph_oversampling/100));
        % params.gen.fov(3) =  round(params.gen.fov(3).*(1+(params.gen.ph_oversampling/100)) / params.gen.res(3)/2)*params.gen.res(3)*2;
        % params.gen.del_k(3) = (1/params.gen.fov(3))*params.gen.kz; 
        % params.gen.n(3) = round(params.gen.fov(3)/params.gen.res(3)/params.gen.kz/2)*2;
        % Just increasing del_k and n..
        params.gen.del_k(3) = params.gen.del_k(3)/(1+(params.gen.ph_oversampling/100)); 
        params.gen.n(3) = round(round(params.gen.n(3)*(1+(params.gen.ph_oversampling/100))/params.gen.kz)/2)*2;
    else
        params.gen.n(3) = round(round(params.gen.n(3)/params.gen.kz)/2)*2; 
    end
    % Fat sat angle
    if params.gen.fs_angle == 0
        params.gen.fs_angle = 110; % 110/80
    end
    
    % Trying to make n multiple of 4,update res. % ToDo: Do I need this?
%     tmp = mod(params.gen.n,4);
%     params.gen.n(1:2) = params.gen.n(1:2)+tmp(1:2);
    % ToDo: Check... Not sure if I do need this...
    % params.gen.n(3) = round(params.gen.n(3)/params.gen.kz/params.gen.pf/2)*2;
    % % Making sure fovz/rz are integers and even
    % if round((params.gen.fov(3)/params.gen.kz)*1000,9)/round((params.gen.fov(3)/params.gen.kz)*1000) ~= 1
    %     params.gen.fov(3) = round(params.gen.fov(3)*1000/2/params.gen.kz)*2*params.gen.kz*1e-3;
    % end
    % if params.gen.seq == 3; params.gen.ro_type = 'c'; end   % if fieldmap, cartesian
    if params.gen.seq == 2; params.gen.ro_type = 's'; end   % if ABC, Spiral
    if params.mt.mt == 0; params.mt.bold = 0; end           % if no MT, no ref BOLD  
    if params.gen.ro_type == 'c'
        params.gen.seg = params.epi.seg;
    elseif params.gen.ro_type == 's'
        params.gen.seg = params.spi.interl;
    end

end
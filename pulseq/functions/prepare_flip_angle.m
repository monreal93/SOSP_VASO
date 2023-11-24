function params = prepare_flip_angle(gx, params)

    % Flip angle 
    % tr_tmp Rough estimate of TR, to calculate Ernst Angle
    if params.gen.ro_type == 'c'
        n_lines = ceil(params.gen.n(2)/params.epi.ry)-((round(params.gen.n(2)/params.epi.ry)-round(params.gen.n(2)/params.epi.ry*params.epi.pf)));
        tr_tmp = (mr.calcDuration(gx)*(n_lines+3))+mr.calcDuration(gx_spoil)+params.gen.te; % +3 of navigators
    elseif params.gen.ro_type == 's'
        tr_tmp = (mr.calcDuration(gx)*params.gen.echos)+2.6e-3+params.gen.te;
    end
    
    if params.gen.fa == 0
        % Ernst angle
        params.gen.fa(1) = acos(exp(-(tr_tmp)/(params.gen.ernst_t1)));
    else
        % Custom angle
        params.gen.fa(1) = params.gen.fa*pi/180;
    end


end
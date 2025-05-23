function [gx_spoil,gy_spoil,gz_spoil] = prepare_spoilers(params,grad_area)
    lims = params.gen.lims;

    max_grad_spoil = mr.convert(params.gen.max_grad_spoil,'mT/m','Hz/m');
    max_sr_spoil = mr.convert(params.gen.max_sr_spoil,'mT/m/ms','Hz/m/s');

    % Set max spoiler gradient and SR
    if params.gen.max_grad_spoil == 0
        params.gen.max_grad_spoil = lims.maxGrad;
    end
    if params.gen.max_sr_spoil == 0
        params.gen.max_sr_spoil = lims.maxSlew;
    end

    gx_spoil=mr.makeTrapezoid('x',lims,'Area',-grad_area*1.5,... 
        'MaxGrad',max_grad_spoil,'MaxSlew',max_sr_spoil);
    gy_spoil=mr.makeTrapezoid('y',lims,'Area',grad_area*1.5,...
        'MaxGrad',max_grad_spoil,'MaxSlew',max_sr_spoil);
    gz_spoil=mr.makeTrapezoid('z',lims,'Area',-grad_area*1.5,...
        'MaxGrad',max_grad_spoil,'MaxSlew',max_sr_spoil);

end
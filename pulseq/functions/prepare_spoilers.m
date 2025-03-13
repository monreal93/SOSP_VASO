function [gx_spoil,gy_spoil,gz_spoil] = prepare_spoilers(params,grad_area)
    lims = params.gen.lims;

%     gx_spoil=mr.makeTrapezoid('x','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);
%     gy_spoil=mr.makeTrapezoid('y','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);
%     gz_spoil=mr.makeTrapezoid('z','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);

    gx_spoil=mr.makeTrapezoid('x',lims,'Area',-grad_area*1.5,'Duration',2.6e-3);
    gy_spoil=mr.makeTrapezoid('y',lims,'Area',grad_area*1.5,'Duration',2.6e-3);
    gz_spoil=mr.makeTrapezoid('z',lims,'Area',-grad_area*1.5,'Duration',2.6e-3);

end
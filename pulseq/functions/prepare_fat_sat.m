function [rf_fs,gx_fs,gy_fs,gz_fs] = prepare_fat_sat(params, sat_ppm)
    lims = params.gen.lims;

    sat_freq = sat_ppm*1e-6*lims.B0*lims.gamma;
    
    % What is the duration I want for the RF pulses??? (2.56 or 5ms)
    rf_fs = mr.makeGaussPulse(params.gen.fs_angle*pi/180,'system',lims,'Duration',2.56e-3,...
            'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
    
    % AMM: ToDo: Check correct values for spoiling gradient
    grad_moment = 25e-6; % T/m * s   (Got this moment from POET simulations of SIEMES EPI)
    grad_ampl = 8e-3;   % T/m       (SIEMENS use 8mT/m)             
    grad_dur = grad_moment/grad_ampl; % s
    grad_dur = round(grad_dur/lims.gradRasterTime)*lims.gradRasterTime;
    grad_area = mr.convert(grad_moment,'T/m/s','Hz/m/s');
    
%     % Using specific amplitude/duration
%     gx_fs = mr.makeTrapezoid('x',lims,'MaxGrad',grad_ampl*lims.gamma,'Area',grad_area,'Duration',grad_dur);
%     gy_fs = mr.makeTrapezoid('y',lims,'MaxGrad',grad_ampl*lims.gamma,'Area',-grad_area,'Duration',grad_dur);
%     gz_fs = mr.makeTrapezoid('z',lims,'MaxGrad',grad_ampl*lims.gamma,'Area',grad_area,'Duration',grad_dur);

    % Using system limits
    gx_fs = mr.makeTrapezoid('x',lims,'Area',grad_area,'Duration',grad_dur);  % Original (1e-3 dur)
    gy_fs = mr.makeTrapezoid('y',lims,'Area',-grad_area,'Duration',grad_dur);
    gz_fs = mr.makeTrapezoid('z',lims,'Area',grad_area,'Duration',grad_dur);

%     gx_fs.delay = mr.calcDuration(rf_fs);
%     gy_fs.delay = mr.calcDuration(rf_fs);
%     gz_fs.delay = mr.calcDuration(rf_fs);

end
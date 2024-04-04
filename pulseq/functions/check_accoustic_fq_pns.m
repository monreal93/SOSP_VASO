function check_accoustic_fq_pns(seq,params, grad_file)
    lims = params.gen.lims;
    % t     time of samples in microseconds, y     sample values, Fs    sampling rate in Hz
    Fs = 1/lims.adcRasterTime/100;
    gradients = seq.gradient_waveforms1()/lims.gamma*1000;   % Full sequence
    gradients = gradients(1:2,:);               % taking only gx and gy

    % Full sequence
    time = linspace(0,seq.duration(),size(gradients,2)).*1e6; % Full sequence
    
    % Only 1 readout
    gradients = gradients(:,1:floor(length(gradients)/params.gen.n_ov(3)));   % 1 readout
    time = time(1:length(gradients));   % 1 readout

    gaxes = ['X' 'Y'];
    for i=1:length(gaxes)   
        gradFreqPlot_pulseq(time,gradients(i,:),Fs,gaxes(i),params.gen.field_strength);
    end
    
    seq.calcPNS(grad_file);

end
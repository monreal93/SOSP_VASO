function check_accoustic_fq_pns(seq,params,seq_t0)
    lims = params.gen.lims;
    % t     time of samples in microseconds, y     sample values, Fs    sampling rate in Hz
    Fs = 1/lims.adcRasterTime/100;

    gradients = seq.gradient_waveforms1()/lims.gamma*1000;   % Full sequence
    tmp = (round(seq_t0)*1e-2 / params.gen.lims.adcRasterTime)+1; % Gradients start point
    gradients = gradients(1:2,tmp:end);               % taking only gx and gy

    % Full sequence (all readouts)
%     time = linspace(0,seq.duration()-seq_t0,size(gradients,2)).*1e6; % Full sequence
    time = (0:size(gradients,2)-1).*params.gen.lims.adcRasterTime.*1e8;
    
    % Only 1 readout
    gradients = gradients(:,1:floor(length(gradients)/params.gen.n_ov(3)));   % 1 readout
    time = time(1:length(gradients));   % 1 readout

    gaxes = ['X' 'Y'];
    for i=1:length(gaxes)   
        [~, ~, ~, ~, bpass(i)] = gradFreqPlot_pulseq(time,gradients(i,:),Fs,gaxes(i),params.gen.field_strength);
    end
    
    % For safety reasons if test fails, the script will give a warning .. 
    if ~bpass(1) || ~bpass(2)
        warndlg(sprintf('########################################## \n Sequence failed resonance test, if you dont feel confident running this sequence, modify Spiral/EPI gradient Amplitde and Slew Rate and try again \n##########################################'))
    end

    if params.gen.pns_check
        if ~isempty(params.gen.grad_file)
            seq.calcPNS(params.gen.grad_file);
        end
    end
    
end
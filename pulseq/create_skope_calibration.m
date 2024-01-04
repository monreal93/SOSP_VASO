% (c) 2022 Skope Magnetic Resonance Technologies AG

clear all
close all
clc

clear all; clc
% Change directory to path where sosp_vaso git has been cloned
cd /home/amonreal/Documents/PhD/tools/skope-pulseq/

% %% Check if Pulseq module has been added
% if not(isfolder('pulseq/matlab'))
%     error("Please run 'git submodule init' and 'git submodule update' to get the latest Pulseq scripts.")
% end

%% Add Pulseq, sequences and methods
% addpath('pulseq/matlab')
addpath('methods')
addpath('sequences')
addpath(genpath("/home/amonreal/Documents/PhD/tools/pulseq/"))

%% Select scanner
scannerType = "Siemens Terra 7T SC72CD";

%% Create off-resonance and position calibration sequence
opc = skope_offresAndPosCalib(scannerType);

% Plot sequence
timeRange = [0 5];
opc.plot(timeRange);

% Test sequence
% opc.test();

namestr = strcat('./pulseq/skope_calibration');
opc.seq.write(strcat(namestr,'.seq'));

save(strcat(namestr,'_params.mat'),'opc')
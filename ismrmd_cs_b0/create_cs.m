%% Adding paths
path_tmp = pwd;
if contains(path_tmp,"amonreal")
% Local paths
    cd /home/amonreal/Documents/PhD/PhD_2022/
    addpath(genpath('/home/amonreal/ismrmrd'))
    addpath(genpath('/home/amonreal/Documents/PhD/tools/matlab/functions/general'))
    addpath(genpath('/home/amonreal/Documents/PhD/tools/matlab/functions/spiral_VASO'))
    addpath(genpath('/home/amonreal/Documents/PhD/tools/bart-0.7.00'))
    addpath(genpath('/home/amonreal/Documents/PhD/tools/mapvbvd-master/'))
    setenv('TOOLBOX_PATH','/home/amonreal/Documents/PhD/tools/bart-0.7.00')
elseif contains(path_tmp,"/mnt/")
    % Beast paths
    cd /mnt/5T4/Alejandro/
    addpath(genpath('/mnt/5T3/Alejandro/tools/ismrmrd'))
    addpath(genpath('/mnt/5T3/Alejandro/tools/matlab/functions/general'))
    addpath(genpath('/mnt/5T3/Alejandro/sosp_vaso/ismrmd_cs_b0/functions'))
    addpath(genpath('/mnt/5T3/Alejandro/tools/bart'))
    addpath(genpath('/mnt/5T3/Alejandro/tools/mapvbvd-master'))
    addpath(genpath("/mnt/5T3/Alejandro/tools/mritools_Linux_3.5.4"))
    setenv('TOOLBOX_PATH','/mnt/5T3/Alejandro/tools/bart')
end

cd ./sosp_vaso

folder = '12062023_sk_abc';
cs_b0_file = 'b0_a01_fieldmap';
scan = 'abc_01_MT_NOfs_NOrfsp';
% repetitions = 1; %4,50        % AMM: ToDo: find a way to get this param from somewhere

% Reading some parameters from Pulseq
load(sprintf('./data/%s/acq/%s_params.mat',folder,scan));

%% Some parameters
params.ncc = 0;              	   % Coil compression coils... 0 for no compression

%% Adding extra parameters from Pulseq
params.slices = params.gen.n_ov(3);
if params.gen.field_strength == 7
    params.ch = 32;  % 7T=32 , 9.4T=31
elseif params.gen.field_strength == 9
    params.ch = 31;  %  9.4T=31, Tubingen coil
%     params.ch = 32;   % 9T visual coil
elseif params.gen.field_strength == 7i
   params.ch = 71;
end
params.nx = params.gen.ro_samples;
params.ny = 1; 
params.rz = params.gen.kz;
params.pf = params.gen.pf;
params.repetitions = repetitions;         % AMM: ToDo: find a way to get this param from somewhere
params.fov = params.gen.fov;  % in mm
params.mtx_s = params.gen.n_ov;
% params.mtx_s(3) = params.mtx_s(3)/params.gen.kz;
params.res = params.gen.res;
params.b0_file = cs_b0_file(4:6);

%% Read b0_map .dat file
[~,twix_params_b0] = fn_read_twix(folder,cs_b0_file,params);
% Getting B0 twix parametersfrom previous save
if exist(sprintf('./data/%s/acq/%s_twix_params_b0.mat',folder,params.b0_file))
    load(sprintf('./data/%s/acq/%s_twix_params_b0.mat',folder,params.b0_file));
else
    fprintf('--- B0 twix parameters havent been read, please run script again and select Yes to read the data \n');
    return
end
params.twix_params_b0 = twix_params_b0;

% %% Read other scans
% [twix_params,~] =  fn_read_twix(folder,scan,params);
% % Getting B0 twix parametersfrom previous save
% if exist(sprintf('./data/%s/acq/%s_twix_params.mat',folder,scan))
%     load(sprintf('./data/%s/acq/%s_twix_params.mat',folder,scan));
% else
%     fprintf('--- SV twix parameters havent been read, please run script again and select Yes to read the data \n');
%     return
% end
% params.twix_params = twix_params;

% %% Create ISMRMD files, one per each repetition/dynamic
% scan_orig = scan;
% if params.gen.echos > 1
%    for i_ech=1:params.gen.echos
%     scan = sprintf('%s_e%i',scan_orig,i_ech);
%     fn_create_ismrmd2(folder,scan,params);
%    end
% else
%     fn_create_ismrmd2(folder,scan,params);
% end
% scan = scan_orig;

%% Generate Coil Sensitivities
fn_coil_sensitivities(folder,scan,cs_b0_file,params);

% %%  Generate B0 map
% fn_get_b0_romeo(folder,scan,cs_b0_file,params);

cd ..

%% Final message
fprintf('--------------------------------------------------------------\n');
fprintf('Sensitivity Maps been created and saved... \n');
fprintf('--------------------------------------------------------------\n');


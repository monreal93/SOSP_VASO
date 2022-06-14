% To do:
% - C1: Check with Skope data if I need to change the sign
% - C2: Check if I need to FFT_1D in partition direction, with data from
% scanner
% - Fix parameter sample_time_us = to get it from file or somewhere else

function fn_create_ismrmd(folder,scan,params)
    
    fprintf('--- Creating ISMRMD files for scan %s ... \n',scan);
    
    % Change parameters if 2D
    if params.is2d == 1
        params.fov = [params.fov(1) params.fov(2) params.fov(3)/params.slices];
        params.slices = 1;
    end

    if params.is2d == 1
        filename_v = sprintf('./data/%s/ismrmd/%s_v_r%i_2d.h5',folder,scan,params.rep_to_save);
        filename_b = sprintf('./data/%s/ismrmd/%s_b_r%i_2d.h5',folder,scan,params.rep_to_save);
    else
        filename_v = sprintf('./data/%s/ismrmd/%s_v_r%i_3d.h5',folder,scan,params.rep_to_save);
        filename_b = sprintf('./data/%s/ismrmd/%s_b_r%i_3d.h5',folder,scan,params.rep_to_save);
    end
    warning('off')
        delete(filename_v);
        delete(filename_b);
    warning('on')
    
    % Create an empty ismrmrd dataset
    if exist(filename_v,'file')
        error(['File ' filename_v ' already exists.  Please remove first'])
    end
    dset_v = ismrmrd.Dataset(filename_v);
    dset_b = ismrmrd.Dataset(filename_b);
    acqData = ismrmrd.Acquisition(1);

    %% Load Data
    % Load Kdata
    load(['./data/' folder '/raw/' scan '_ks_vaso.mat']);
    ks_vaso = squeeze(ks_vaso(:,:,params.rep_to_save,:));
    load(['./data/' folder '/raw/' scan '_ks_bold.mat']);
    ks_bold = squeeze(ks_bold(:,:,params.rep_to_save,:));

    % Load Trajectory
    if params.traj == 3
        load(['./data/' folder '/acq/' scan '_ks_traj_sk.mat']);
        % C1: Not sure about this, need to check with Skope data
        ks_traj.kx = ks_traj.kx.*-1;
        ks_traj.ky = ks_traj.ky.*-1;
        ks_traj.kz = ks_traj.kz.*-1;
    elseif params.traj == 1
        load(['./data/' folder '/acq/' scan '_ks_traj_nom.mat']);
        ks_traj.ky = ks_traj.ky.*-1;        % Swaping nominal trajectory to match scaner one
    elseif params.traj == 2
        load(['./data/' folder '/acq/' scan '_ks_traj_nom_poet.mat']);
        ks_traj.kx = ks_traj.kx.*-1;        % Swaping nominal trajectory to match scaner one
        ks_traj.ky = ks_traj.ky.*-1;        % Swaping nominal trajectory to match scaner one
    end

%     % Cropping trajctory to center and edges of kspace, if its from Skope is
%     % should be cropped manually before
%     if params.traj ~= 3
%         % very ugly way to find the start of the spiral
%         tmp = diff(diff(ks_traj.kx(1,1,:)));
%         tmp = tmp<1e-7;
%         tmp = tmp == 1; tmp = single(tmp);
%         tmp(tmp==0) = 2; tmp(tmp==1) = 0;
%         tmp = find(tmp); 
%         idx = min(tmp)-1;
%         ks_traj.kx = ks_traj.kx(:,:,idx:end-1);
%         ks_traj.ky = ks_traj.ky(:,:,idx:end-1);
%         ks_traj.kz = ks_traj.kz(:,:,idx:end-1);   
%     end

    % Normalize and interpolate trajectory
    ks_traj = norm_interp_traj(ks_traj.kx,ks_traj.ky,ks_traj.kz,params.nx*params.spi.interl);

%     %%%%%%%%% Temp, Trying to match RECON as in gpuNUFFT
%     % shifting BOLD slices
%     tmp = FFT_3D(ks_bold,'image');
%     tmp = circshift(tmp,[0,8,0,0]);
%     tmp = FFT_3D(tmp,'kspace');
% %     tmp = FFT_1D(tmp,'image',2);  % Putting it in the hyb dom to work as 2D 
%     ks_bold = tmp;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % FFT shift data
    % ks_vaso = fft_shift_2D(ks_vaso,ks_traj.kx,0,ks_traj.ky,round(params.twix_params_b0.shift/params.res(2)./1000)); % 14,31
    % ks_bold = fft_shift_2D(ks_bold,ks_traj.kx,0,ks_traj.ky,round(params.twix_params_b0.shift/params.res(2)./1000));

    % Subseting data if is 2d
    if params.is2d==1 
        if params.traj== 1
            ks_traj.kx = ks_traj.kx(:,params.slice_to_save);
            ks_traj.ky = ks_traj.ky(:,params.slice_to_save);
            ks_traj.kz = ks_traj.kz(:,params.slice_to_save);
        end
        % AMM: Here I don't know why for sv_1 ks_vaso, I need to do FFT..
        ks_vaso = FFT_1D(ks_vaso,'image',2);
        ks_vaso = ks_vaso(:,params.slice_to_save,:,:);
        ks_bold = FFT_1D(ks_bold,'image',2);
        ks_bold = ks_bold(:,params.slice_to_save,:,:);
    end

    % Permuting to make diemnsions compatible with ISMRM reader in Julia
    ks_traj.kx = permute(ks_traj.kx,[3,2,1]);
    ks_traj.ky = permute(ks_traj.ky,[3,2,1]);
    ks_traj.kz = permute(ks_traj.kz,[3,2,1]);

    % %%%%%%%%% Temp, Trying to match RECON as in gpuNUFFT
    % % shifting BOLD slices
    % tmp = FFT_3D(ks_bold,'image');
    % tmp = circshift(tmp,[0,8,0,0]);
    % tmp = FFT_3D(tmp,'kspace');
    % tmp = FFT_1D(tmp,'image',2);  % Putting it in the hyb dom to work as 2D 
    % ks_bold = tmp;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % % Subseting data if is 2d
    % if params.is2d==1
    %     ks_traj.kx = ks_traj.kx(:,slice_to_save);
    %     ks_traj.ky = ks_traj.ky(:,slice_to_save);
    %     ks_traj.kz = ks_traj.kz(:,slice_to_save);
    % %     ks_vaso = FFT_1D(ks_vaso,'image',2);
    %     ks_vaso = ks_vaso(:,slice_to_save,:);
    % %     ks_bold = FFT_1D(ks_bold,'image',2);
    %     ks_bold = ks_bold(:,slice_to_save,:);
    % end

    %% Setting header parameters
    % Set the header elements that don't change
    acqData.head.version(:) = 1;
    acqData.head.number_of_samples(:) = params.nx*params.spi.interl;
    acqData.head.center_sample(:) = floor(params.nx*params.spi.interl/2);
    acqData.head.active_channels(:) = params.ch;
    acqData.head.read_dir  = repmat([1 0 0]',[1 params.ny]);
    acqData.head.phase_dir = repmat([0 1 0]',[1 params.ny]);
    acqData.head.slice_dir = repmat([0 0 1]',[1 params.ny]);
    acqData.head.sample_time_us = (1.535/2);  % In seconds (should be in seconds for MRIReco)

    % AMM:
    if params.is2d == 1
        acqData.head.trajectory_dimensions = 2;
    else
        acqData.head.trajectory_dimensions = 3;
    end

    %% Adding new data
    ks_vaso = double(ks_vaso);
    ks_bold = double(ks_bold);

    % VASO
    for i=1:floor(params.slices/params.rz*params.pf)
        if params.is2d == 1
            traj = [squeeze(ks_traj.kx(1,i,:)) squeeze(ks_traj.ky(1,i,:))];
        else
            traj = [squeeze(ks_traj.kx(1,i,:)) squeeze(ks_traj.ky(1,i,:)) squeeze(ks_traj.kz(1,i,:))];
        end
        data = squeeze(ks_vaso(:,i,:));

        if i == 1
            tmp_data = data;
            tmp_traj = traj';
        else
            tmp_data = cat(1,tmp_data,data);
            tmp_traj = cat(2,tmp_traj,traj');
        end
    %     acqData.data{i} = data;
    %     acqData.traj{i} = traj';
    end
    acqData.data{1} = tmp_data;
    acqData.traj{1} = tmp_traj;
    dset_v.appendAcquisition(acqData);

    % BOLD
    for i=1:floor(params.slices/params.rz*params.pf)
        if params.is2d == 1
            traj = [squeeze(ks_traj.kx(1,i,:)) squeeze(ks_traj.ky(1,i,:))];
        else
            traj = [squeeze(ks_traj.kx(1,i,:)) squeeze(ks_traj.ky(1,i,:)) squeeze(ks_traj.kz(1,i,:))];
        end
        data = squeeze(ks_bold(:,i,:));

        if i == 1
            tmp_data = data;
            tmp_traj = traj';
        else
            tmp_data = cat(1,tmp_data,data);
            tmp_traj = cat(2,tmp_traj,traj');
        end
    %     acqData.data{i} = data;
    %     acqData.traj{i} = traj';
    end
    acqData.data{1} = tmp_data;
    acqData.traj{1} = tmp_traj;
    dset_b.appendAcquisition(acqData);

    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Fill the xml header %
    %%%%%%%%%%%%%%%%%%%%%%%%
    % We create a matlab struct and then serialize it to xml.
    % Look at the xml schema to see what the field names should be

    header = [];

    % Experimental Conditions (Required)
    header.experimentalConditions.H1resonanceFrequency_Hz = 298000000; % 7T

    % Acquisition System Information (Optional)
    header.acquisitionSystemInformation.systemVendor = 'ISMRMRD Labs';
    header.acquisitionSystemInformation.systemModel = 'Virtual Scanner';
    header.acquisitionSystemInformation.receiverChannels = params.ch;

    % The Encoding (Required)
    header.encoding.trajectory = 'spiral';
    header.encoding.encodedSpace.fieldOfView_mm.x = params.fov(1);
    header.encoding.encodedSpace.fieldOfView_mm.y = params.fov(2);
    header.encoding.encodedSpace.fieldOfView_mm.z = params.fov(3);
    header.encoding.encodedSpace.matrixSize.x = size(traj,1);
    header.encoding.encodedSpace.matrixSize.y = 1;
    header.encoding.encodedSpace.matrixSize.z = params.slices;
    % Recon Space
    % (in this case same as encoding space)
    header.encoding.reconSpace = header.encoding.encodedSpace;
    % Encoding Limits
    header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_0.maximum = size(traj,1)-1;
    header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(size(traj,1)/2);
    header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_1.maximum = size(traj,2)-1;
    header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(size(traj,3)/2);

    header.encoding.encodingLimits.repetition.minimum = 0;
    header.encoding.encodingLimits.repetition.maximum = params.repetitions-1;
    header.encoding.encodingLimits.repetition.center = 0;

    header.encoding.encodingLimits.slice.minimum = 0;
    header.encoding.encodingLimits.slice.maximum = params.slices;
    header.encoding.encodingLimits.slice.center = floor(params.slices/2);
    % header.encoding.encodingLimits.repetition.maximum = 1;

    %% Serialize and write to the data set
    xmlstring = ismrmrd.xml.serialize(header);
    dset_v.writexml(xmlstring);
    dset_b.writexml(xmlstring);

    dset_v.close();
    dset_b.close();
    
    fprintf('ISMRMD file repetition %i saved in path %s/ismrmd/ \n',params.rep_to_save, folder);

        if params.plot
            % Plotting some stuff
            % 3D traj
            tmp = rssq(ks_bold,3);
            tmp = tmp(:).*1e6;
            if params.is2d == 0
                figure;scatter3(ks_traj.kx(:),ks_traj.ky(:),ks_traj.kz(:),[],tmp,'.')
            end

            % 2D traj
            tmp = rssq(ks_bold(:,1,:),3);
            tmp = tmp(:).*1e6;
            figure;scatter(ks_traj.kx(:,1),ks_traj.ky(:,1),[],tmp,'.'); title('K-space/trajectory sample')
        end
              
end

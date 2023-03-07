% To do:
% This function is to be called inside fn_read_twix, it won't read the
% raw data from a .mat file, it will recieve k_vaso and k_bold as
% parameters

function fn_create_ismrmd5(folder,scan,params,ks_vaso_all,ks_bold_all)
    % This function creates one mrd file per repetition and per slice...
    fprintf('--- Creating ISMRMD files for scan %s ... \n',scan);
    
    % Change parameters if 2D
    if params.is2d == 1
        params.fov = [params.fov(1) params.fov(2) params.fov(3)/params.slices];
        % params.slices = 1;
    end
    
    % Check trajectory type
    if params.traj == 1
        traj_name = 'nom';
    elseif params.traj == 2
        traj_name = 'poet';
    elseif params.traj == 3
        traj_name = 'sk';
    end
    
    % Checking if files exist, it checks for the first one only
    if contains(scan,'sv') ||  contains(scan,'cv')
        if params.is2d == 1
            filename = sprintf('./data/%s/ismrmd/2d/%s_v_r1_sl1_2d_%s.h5',folder,scan,traj_name);
        else
            filename = sprintf('./data/%s/ismrmd/3d/%s_v_r1_3d_%s.h5',folder,scan,traj_name);
        end
    elseif contains(scan,'abc')
        if params.is2d == 1
            filename = sprintf('./data/%s/ismrmd/2d/%s_r1_sl1_2d_%s.h5',folder,scan,traj_name);
        else
            filename = sprintf('./data/%s/ismrmd/3d/%s_r1_3d_%s.h5',folder,scan,traj_name);
        end
    end
    tmp = 'y';
    if exist(filename) > 1
        prompt = 'ISMRMD files already exist, do you want to replace it: type y/n [y]: ';
        tmp = input(prompt,'s');
    end
    
    % Gradient delay
    % dw = params.gen.t_vector(2)-params.gen.t_vector(1);
    % g_delay = params.twix_params.grad_delay./dw*1e-6;
    
    % Check selection:
    if tmp == 'n'
            fprintf('B0 map not generated \n');  
            return
    elseif tmp == 'y'
        % Load K-space data
%         if contains(scan,'sv')  ||  contains(scan,'cv')
%             load(['./data/' folder '/raw/' scan '_ks_vaso.mat']);
%             ks_vaso_all =ks_vaso;
%             load(['./data/' folder '/raw/' scan '_ks_bold.mat']);
%             ks_bold_all = ks_bold;
%         elseif contains(scan,'abc')
%             load(['./data/' folder '/raw/' scan '_ks_abc.mat']);
%             ks_abc_all = ks_abc;
%         end
        
        % If is 2D, repeat this for every slice
        if params.is2d
            partitions = params.slices;
            params.slice_to_save = 1:params.slices;
            last_partition = 1;
        else
            last_partition = floor(params.slices/params.rz*params.pf);
            partitions = 1;
        end
        
        % Loop over repetitions
        for j = 1:params.repetitions
            for k=1:partitions % slices
                % Loop over slices if 2D
                params.rep_to_save = j;
                warning('off')
                if contains(scan,'sv')  ||  contains(scan,'cv')
                    if params.is2d == 1
                        filename_v = sprintf('./data/%s/ismrmd/2d/%s_v_r%i_sl%i_2d_%s.h5',folder,scan,j,k,traj_name);
                        filename_b = sprintf('./data/%s/ismrmd/2d/%s_b_r%i_sl%i_2d_%s.h5',folder,scan,j,k,traj_name);
                    else
                        filename_v = sprintf('./data/%s/ismrmd/3d/%s_v_r%i_3d_%s.h5',folder,scan,params.rep_to_save,traj_name);
                        filename_b = sprintf('./data/%s/ismrmd/3d/%s_b_r%i_3d_%s.h5',folder,scan,params.rep_to_save,traj_name);
                    end
                    delete(filename_v);
                    delete(filename_b);
                elseif contains(scan,'abc')
                    if params.is2d == 1
                        filename = sprintf('./data/%s/ismrmd/2d/%s_r%i_sl%i_2d_%s.h5',folder,scan,j,k,traj_name);
                    else
                        filename = sprintf('./data/%s/ismrmd/3d/%s_r%i_3d_%s.h5',folder,scan,params.rep_to_save,traj_name);
                    end
                    delete(filename);
                end
                warning('on')

            %     % Create an empty ismrmrd dataset
            %     if exist(filename_v,'file')
            %         error(['File ' filename_v ' already exists.  Please remove first'])
            %     end
            %     dset_v = ismrmrd.Dataset(filename_v);
            %     dset_b = ismrmrd.Dataset(filename_b);
                acqData = ismrmrd.Acquisition(1);

                if contains(scan,'sv')  ||  contains(scan,'cv')
                    dset_v = ismrmrd.Dataset(filename_v);
                    dset_b = ismrmrd.Dataset(filename_b);
                    ks_vaso = squeeze(ks_vaso_all(:,:,params.rep_to_save,:));
                    ks_bold = squeeze(ks_bold_all(:,:,params.rep_to_save,:));
                elseif contains(scan,'abc')
                    dset_abc = ismrmrd.Dataset(filename);
                    ks_abc = squeeze(ks_abc_all(:,:,params.rep_to_save,:));
                end

                % Load Trajectory
                if params.traj == 3
                    load(['./data/' folder '/acq/' scan '_ks_traj_sk.mat']);
                    % C1: Not sure about this, need to check with Skope data
                    % ks_traj.kx = ks_traj.kx.*-1;
                    ks_traj.ky = ks_traj.ky.*-1;
                    % ks_traj.kz = ks_traj.kz.*-1;
                elseif params.traj == 1
                    load(['./data/' folder '/acq/' scan '_ks_traj_nom.mat']);
                    ks_traj.ky = ks_traj.ky.*-1;        % Swaping nominal trajectory to match scaner one
                elseif params.traj == 2
                    load(['./data/' folder '/acq/' scan '_ks_traj_poet.mat']);
%                     ks_traj.kx = ks_traj.kx.*-1;        % Swaping nominal trajectory to match scaner one
                    ks_traj.ky = ks_traj.ky.*-1;        % Swaping nominal trajectory to match scaner one
                end

                %%% Temp: Discaring corrupted samples
                st_crop = find(ks_traj.kx(:,1));
                st_crop = st_crop(1); st_crop = floor(st_crop/10)*10+1;
                ks_traj.kx = ks_traj.kx(st_crop:end,:);
                ks_traj.ky = ks_traj.ky(st_crop:end,:);
                ks_traj.kz = ks_traj.kz(st_crop:end,:);
                if params.gen.seq == 1
                    ks_vaso = ks_vaso(st_crop:end,:,:);
                    ks_bold = ks_bold(st_crop:end,:,:);
                elseif params.gen.seq == 2
                    ks_abc = ks_abc(st_crop:end,:,:);
                end
                % Updating params...
                params.nx = length(ks_traj.kx);
                params.gen.t_vector = params.gen.t_vector(st_crop:end);
                params.spi.ro_samples = length(ks_traj.kx);
                save(sprintf('./data/%s/acq/%s_params.mat',folder,scan),'params');
                %%%%%
                
                
%                 % Correcting for gradient Delay
%                 g_delay = 89e-9;
%                 kx = interp1(params.gen.t_vector,ks_traj.kx,params.gen.t_vector+g_delay,'linear','extrap');
%                 ky = interp1(params.gen.t_vector,ks_traj.ky,params.gen.t_vector+g_delay,'linear','extrap');
%                 kz = interp1(params.gen.t_vector,ks_traj.kz,params.gen.t_vector+g_delay,'linear','extrap');
%                 ks_traj.kx = kx;
%                 ks_traj.ky = ky;
%                 ks_traj.kz = kz;

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

                % Temp: FFT shift data, for now only spiral, not sure if
                % needed for cartesian
                if params.gen.ro_type == 's'
                    if params.gen.seq == 1
                        ks_vaso = fft_shift_2D(ks_vaso,ks_traj.kx,round(params.twix_params_b0.shift/params.res(2)./10000),ks_traj.ky,0); % 14,31
                        ks_bold = fft_shift_2D(ks_bold,ks_traj.kx,round(params.twix_params_b0.shift/params.res(2)./10000),ks_traj.ky,0);
                    elseif params.gen.seq == 2
                        ks_abc = fft_shift_2D(ks_abc,ks_traj.kx,round(params.twix_params_b0.shift/params.res(2)./10000),ks_traj.ky,0);
                    end
                end
                
                % Subseting data if is 2d
                if params.is2d==1 
                    if params.traj== 1
                        ks_traj.kx = ks_traj.kx(:,params.slice_to_save(k));
                        ks_traj.ky = ks_traj.ky(:,params.slice_to_save(k));
                        ks_traj.kz = ks_traj.kz(:,params.slice_to_save(k));
                    end
                    if contains(scan,'sv')  ||  contains(scan,'cv')
                        % AMM: Here I don't know why for sv_1 ks_vaso, I need to do FFT..
                        ks_vaso = FFT_1D(ks_vaso,'image',2);
                        ks_vaso = ks_vaso(:,params.slice_to_save(k),:,:);
                        ks_bold = FFT_1D(ks_bold,'image',2);
                        ks_bold = ks_bold(:,params.slice_to_save(k),:,:);
                    elseif contains(scan,'abc')
                        ks_abc = FFT_1D(ks_abc,'image',2);
                        ks_abc = ks_abc(:,params.slice_to_save(k),:,:);
                    end
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
                if contains(scan,'sv')  ||  contains(scan,'cv')
                    ks_vaso = double(ks_vaso);
                    ks_bold = double(ks_bold);

                    % VASO
                    for i=1:last_partition
                        if params.is2d == 1
                            traj = [squeeze(ks_traj.kx(1,1,:)) squeeze(ks_traj.ky(1,1,:))];
                            data = squeeze(ks_vaso(:,1,:));
                        else
                            traj = [squeeze(ks_traj.kx(1,i,:)) squeeze(ks_traj.ky(1,i,:)) squeeze(ks_traj.kz(1,i,:))];
                            data = squeeze(ks_vaso(:,i,:));
                        end
                        

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
                    for i=1:last_partition
                        if params.is2d == 1
                            traj = [squeeze(ks_traj.kx(1,1,:)) squeeze(ks_traj.ky(1,1,:))];
                            data = squeeze(ks_bold(:,1,:));
                        else
                            traj = [squeeze(ks_traj.kx(1,i,:)) squeeze(ks_traj.ky(1,i,:)) squeeze(ks_traj.kz(1,i,:))];
                            data = squeeze(ks_bold(:,i,:));
                        end

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
                elseif contains(scan,'abc')
                    ks_abc = double(ks_abc);

                    % ABC
                    for i=1:last_partition
                        if params.is2d == 1
                            traj = [squeeze(ks_traj.kx(1,1,:)) squeeze(ks_traj.ky(1,1,:))];
                            data = squeeze(ks_abc(:,1,:));
                        else
                            traj = [squeeze(ks_traj.kx(1,i,:)) squeeze(ks_traj.ky(1,i,:)) squeeze(ks_traj.kz(1,i,:))];
                            data = squeeze(ks_abc(:,i,:));
                        end
                        

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
                    dset_abc.appendAcquisition(acqData);

                end

                %%%%%%%%%%%%%%%%%%%%%%%%
                %% Fill the xml header %
                %%%%%%%%%%%%%%%%%%%%%%%%
                % We create a matlab struct and then serialize it to xml.
                % Look at the xml schema to see what the field names should be

                header = [];

                % Experimental Conditions (Required)
                header.experimentalConditions.H1resonanceFrequency_Hz = 298000000; % 7T

                % Acquisition System Information (Optional)
                header.acquisitionSystemInformation.systemVendor = 'Scannexus';
            %     header.acquisitionSystemInformation.systemModel = 'Virtual Scanner';
                header.acquisitionSystemInformation.receiverChannels = params.ch;

                % The Encoding (Required)
                header.encoding.trajectory = 'spiral';
                header.encoding.encodedSpace.fieldOfView_mm.x = params.fov(1);
                header.encoding.encodedSpace.fieldOfView_mm.y = params.fov(2);
                header.encoding.encodedSpace.fieldOfView_mm.z = params.fov(3);
                header.encoding.encodedSpace.matrixSize.x = size(traj,1);
                header.encoding.encodedSpace.matrixSize.y = 1;
                header.encoding.encodedSpace.matrixSize.z = last_partition;
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
                header.encoding.encodingLimits.slice.maximum = last_partition;
                header.encoding.encodingLimits.slice.center = floor(last_partition/2);
                % header.encoding.encodingLimits.repetition.maximum = 1;

                %% Serialize and write to the data set
                xmlstring = ismrmrd.xml.serialize(header);

                if contains(scan,'sv')  ||  contains(scan,'cv')
                    dset_v.writexml(xmlstring);
                    dset_b.writexml(xmlstring);

                    dset_v.close();
                    dset_b.close();
                elseif contains(scan,'abc')
                    dset_abc.writexml(xmlstring);
                    dset_abc.close();  
                end

                fprintf('ISMRMD file repetition %i, slice %i saved in path %s/ismrmd/ \n',params.rep_to_save,k, folder);

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
        end
    end
end

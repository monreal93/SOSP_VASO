% To do:
% - C1: Check with Skope data if I need to change the sign
% - C2: Check if I need to FFT_1D in partition direction, with data from
% scanner
% - Fix parameter sample_time_us = to get it from file or somewhere else

function fn_create_ismrmd2(folder,scan,params)
  % This function creates ismrmd files including 
  % partitions and one file per repetition
  
    fprintf('--- Creating ISMRMD files for scan %s ... \n',scan);
    
    % Change parameters if 2D
    if params.is2d == 1
%         params.fov = [params.fov(1) params.fov(2) params.fov(3)/params.slices];
        params.fov = [params.fov(1) params.fov(2) params.fov(3)];
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
    
%     % Check DORK
%     if params.gen.dork == 1             % partial DORK
%         dork = '_pDORK';
%     elseif params.gen.dork == 2         % Full DORK
%         dork = '_fDORK';
%     else
%         dork = '';
%     end
    if params.part_dork == 1
        p_dork = '_pDORK';
    else
        p_dork = '';
    end
    if params.rep_dork == 1
        r_dork = '_rDORK';
    else
        r_dork = '';
    end
    
    % Checking if files exist, it checks for the first one only
    if contains(scan,'sv') ||  contains(scan,'cv')
        if params.is2d == 1
            filename = sprintf('./data/%s/ismrmd/2d/%s_v_r1_sl1_2d_%s%s%s.h5',folder,scan,traj_name,p_dork,r_dork);
        else
            filename = sprintf('./data/%s/ismrmd/3d/%s_v_3d_r1_%s%s%s.h5',folder,scan,traj_name,p_dork,r_dork);
        end
    elseif contains(scan,'abc')
        if params.is2d == 1
            filename = sprintf('./data/%s/ismrmd/2d/%s_r1_sl1_2d_%s%s%s.h5',folder,scan,traj_name,p_dork,r_dork);
        else
            filename = sprintf('./data/%s/ismrmd/3d/%s_3d_r1_%s%s%s.h5',folder,scan,traj_name,p_dork,r_dork);
        end
    elseif contains(scan,'sb')
        if params.is2d == 1
            filename = sprintf('./data/%s/ismrmd/2d/%s_r1_sl1_2d_%s%s%s.h5',folder,scan,traj_name,p_dork,r_dork);
        else
            filename = sprintf('./data/%s/ismrmd/3d/%s_3d_r1_%s%s%s.h5',folder,scan,traj_name,p_dork,r_dork);
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
        if contains(scan,'sv')  ||  contains(scan,'cv')
            load(['./data/' folder '/raw/' scan '_ks_vaso.mat']);
            load(['./data/' folder '/raw/' scan '_ks_bold.mat']);
        elseif contains(scan,'abc')
            load(['./data/' folder '/raw/' scan '_ks_abc.mat']);
        elseif contains(scan,'sb')
            load(['./data/' folder '/raw/' scan '_ks_sb.mat']);
        end
        
%         last_partition = floor(params.slices/params.rz*params.pf);
        last_partition = params.mtx_s(3);
        % If is 2D, do FFT in partition direction
        if params.is2d
            if contains(scan,'sv')  ||  contains(scan,'cv')
                ks_vaso = FFT_1D(ks_vaso,'image',2);
                ks_bold = FFT_1D(ks_bold,'image',2);
            elseif contains(scan,'abc')
                ks_abc = FFT_1D(ks_abc,'image',2);
            elseif contains(scan,'sb')
                ks_abc = FFT_1D(ks_sb,'image',2);
            end
        end
        
        % Checking if dataset exists and delete it
        warning('off')
        if contains(scan,'sv')  ||  contains(scan,'cv')
            if params.is2d == 1
                filename_v = sprintf('./data/%s/ismrmd/2d/%s_v_2d_%s.h5',folder,scan,traj_name);
                filename_b = sprintf('./data/%s/ismrmd/2d/%s_b_2d_%s.h5',folder,scan,traj_name);
            else
                filename_v = sprintf('./data/%s/ismrmd/3d/%s_v_3d_%s.h5',folder,scan,traj_name);
                filename_b = sprintf('./data/%s/ismrmd/3d/%s_b_3d_%s.h5',folder,scan,traj_name);
            end
            delete(filename_v);
            delete(filename_b);
        elseif contains(scan,'abc')
            if params.is2d == 1
                filename = sprintf('./data/%s/ismrmd/2d/%s_2d_%s.h5',folder,scan,traj_name);
            else
                filename = sprintf('./data/%s/ismrmd/3d/%s_3d_%s.h5',folder,scan,traj_name);
            end
            delete(filename);
        elseif contains(scan,'sb')
            if params.is2d == 1
                filename = sprintf('./data/%s/ismrmd/2d/%s_2d_%s.h5',folder,scan,traj_name);
            else
                filename = sprintf('./data/%s/ismrmd/3d/%s_3d_%s.h5',folder,scan,traj_name);
            end
            delete(filename);
        end
        warning('on')
        
        % Create data set and Acquisition
        if contains(scan,'sv')  ||  contains(scan,'cv')
            dset_v = ismrmrd.Dataset(filename_v);
            dset_b = ismrmrd.Dataset(filename_b);
        elseif contains(scan,'abc')
            dset_abc = ismrmrd.Dataset(filename);
        elseif contains(scan,'sb')
            dset_sb = ismrmrd.Dataset(filename);
        end
        
        acqData = ismrmrd.Acquisition(params.slices);

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
%             ks_traj.kx = ks_traj.kx.*-1;        % Swaping nominal trajectory to match scaner one
            ks_traj.ky = ks_traj.ky.*-1;        % Swaping nominal trajectory to match scaner one
        end
        
%         %%% Temp: Discaring corrupted samples
%         st_crop = find(ks_traj.kx(:,1));
%         st_crop = st_crop(1); st_crop = floor(st_crop/10)*10+1;
%         ks_traj.kx = ks_traj.kx(st_crop:end,:);
%         ks_traj.ky = ks_traj.ky(st_crop:end,:);
%         ks_traj.kz = ks_traj.kz(st_crop:end,:);
%         if params.gen.seq == 1
%             ks_vaso = ks_vaso(st_crop:end,:,:,:);
%             ks_bold = ks_bold(st_crop:end,:,:,:);
%         elseif params.gen.seq == 2
%             ks_abc = ks_abc(st_crop:end,:,:,:);
%         end
%         % Updating params...
%         params.nx = length(ks_traj.kx);
%         params.gen.t_vector = params.gen.t_vector(st_crop:end);
%         params.gen.ro_samples = length(ks_traj.kx);
%         save(sprintf('./data/%s/acq/%s_params.mat',folder,scan),'params');
%         %%%%%
        
        % Normalize and interpolate trajectory
        if params.traj == 3
            k0 = ks_traj.k0;
        end
        ks_traj = norm_interp_traj(ks_traj.kx,ks_traj.ky,ks_traj.kz,params.nx*params.spi.interl);
                
        % Temp: FFT shift data, for now only spiral, not sure if
        % needed for cartesian
        if params.gen.ro_type == 's'
            fft_shift(1) = params.twix_params.shift(1)./params.gen.res(1)./10000;
            fft_shift(2) = params.twix_params.shift(2)./params.gen.res(2)./10000;
            
            if params.gen.seq == 1
%                 ks_vaso = fft_shift_2D(ks_vaso,ks_traj.kx,round(params.twix_params_b0.shift(1)/params.res(1)./10000),ks_traj.ky,round(params.twix_params_b0.shift(2)/params.res(2)./10000)); % 14,31
%                 ks_bold = fft_shift_2D(ks_bold,ks_traj.kx,round(params.twix_params_b0.shift(1)/params.res(1)./10000),ks_traj.ky,round(params.twix_params_b0.shift(2)/params.res(2)./10000));
                ks_vaso = fft_shift_2D(ks_vaso,ks_traj.kx,fft_shift(1),ks_traj.ky,fft_shift(2)); % 14,31
                ks_bold = fft_shift_2D(ks_bold,ks_traj.kx,fft_shift(1),ks_traj.ky,fft_shift(2));
           
            elseif params.gen.seq == 2
%                 ks_abc = fft_shift_2D(ks_abc,ks_traj.kx,round(params.twix_params_b0.shift/params.res(2)./10000),ks_traj.ky,0);
                ks_abc = fft_shift_2D(ks_abc,ks_traj.kx,fft_shift(1),ks_traj.ky,fft_shift(2));
            end
        end
        
%%%%%%%%%%%%%%%%%%% DEMO: Perform phase correction...
if params.gen.dork ~= 0
    
        % Dirty trick for the moment, to make this part of code work with
        % ABC data...
        if params.gen.seq == 2
            ks_bold = ks_abc;
            ks_vaso = ks_abc;
        end
        
        ks_bold_new = ks_bold;
        ks_vaso_new = ks_vaso;
        % Demodulation of k0
        load('/mnt/5T3/Alejandro/sosp_vaso/data/tmp/k0_sk.mat')
        
        if params.k0_demodulation == 1
            ks_bold_new = ks_bold_new.*exp(1i.*(k0));
            ks_vaso_new = ks_vaso_new.*exp(1i.*(k0));
        end
        
         load('/mnt/5T3/Alejandro/sosp_vaso/data/03032023_sv/tmp/sim_kspace.mat')
    
         % Trying to find the phase difference btw partitions in diff
         % repetitions
         
            
        if params.gen.dork == 2
            load(['./data/' folder '/raw/' scan '_nav.mat']);   
            nav.vaso = squeeze(mean(nav.vaso,4));
            nav.bold = squeeze(mean(nav.bold,4));
        elseif params.gen.dork == 1
            tmp = size(ks_vaso_new);
            nav.vaso = zeros(tmp(1:3));
            nav.bold = zeros(tmp(1:3));
        end

        % Do I want to do separate VASO  and BOLD?
        % Here I try to get the navigator,I look in the kspace traj and
        % take the 100 samples before the spiral is played
        tmp = find(ks_traj.kx(:,1));
        if params.traj == 3
            tmp = find(abs(ks_traj.kx(:,1))>2e-5);
        end
        l_nav = 45; % (240)
        f_nav = 10;%l_nav-99;
        ref_rep = 2;
        ks_tmp_b =  mean(ks_bold,4);
        nav_ri_b = (ks_tmp_b(f_nav:l_nav,(params.gen.n(3)/2)+1,ref_rep)); % Reference repetition (first nav @TE)
        nav_rn_b = nav.bold(:,(params.gen.n(3)/2)+1,ref_rep);             % Reference repetition (last nav @TN)
        ks_tmp_v =  mean(ks_vaso,4);
        nav_ri_v = (ks_tmp_v(f_nav:l_nav,(params.gen.n(3)/2)+1,ref_rep)); % Reference repetition (first nav @TE)
        nav_rn_v = nav.vaso(:,(params.gen.n(3)/2)+1,ref_rep);             % Reference repetition (last nav @TN)
        t = params.gen.t_vector.';
        if params.gen.dork == 2
            tn = t(end)+2e-3;
        else
            tn = 0;
        end
        
        tmp_te  = params.gen.TE;%+2e-3;
    %%%% Partition DORK
         if params.part_dork == 1
             for i = 1:params.repetitions
                %%%%%%% AMM: Temp, trying to compare the phase of the simulated data to the
                 % in-vivo dataa
                 % Read the simulated data
                 load('/mnt/5T3/Alejandro/sosp_vaso/data/03032023_sv/tmp/sim_kspace.mat')
                 % put it in the same space as in-vivo
    %              xx = complex(real(ks).*-1,imag(ks));
    %              xx = xx*exp(1i);
                 ks_nav = ks;
                % Here, I try to correct for the phase difference between
                % partitions, I take as reference the simulated data (xx), it works
                % best with del_phi negative...

                % nav_xx = squeeze(ks_nav(150:250,:,:));
                nav_xx = squeeze(ks_nav(f_nav:l_nav,:,:));
                nav_xx = mean(nav_xx,3);
                nav_xx = angle(nav_xx);

                nav_xx_vaso = squeeze(ks_vaso(f_nav:l_nav,:,i,:));
                nav_xx_vaso = mean(nav_xx_vaso,3);
                nav_xx_vaso = angle(nav_xx_vaso);
                del_phi_xx_v = (nav_xx-nav_xx_vaso);%./params.gen.TE; % Nav ph diff
                del_phi_xx_v = mean(del_phi_xx_v,1);
                del_phi_xx_v = repmat(del_phi_xx_v,[params.gen.ro_samples,1,1,1]);
                del_phi_xx_v = permute(del_phi_xx_v,[1,2,4,3]);
                % del_phi_xx_v = del_phi_xx_v.*-1;
%                 ks_vaso_new_new =  ks_vaso_new.*exp(1i*(del_phi_xx_v));
                ks_vaso_new_new(:,:,i,:) =  ks_vaso_new(:,:,i,:).*exp(1i*(del_phi_xx_v));

                
                nav_xx_bold = squeeze(ks_bold(f_nav:l_nav,:,i,:));
                nav_xx_bold = mean(nav_xx_bold,3);
                nav_xx_bold = angle(nav_xx_bold);
                del_phi_xx_b = (nav_xx-nav_xx_bold);%./params.gen.TE; % Nav ph diff
                del_phi_xx_b = mean(del_phi_xx_b,1);
                del_phi_xx_b = repmat(del_phi_xx_b,[params.gen.ro_samples,1,1,1]);
                del_phi_xx_b = permute(del_phi_xx_b,[1,2,4,3]);
                % del_phi_xx_v = del_phi_xx_v.*-1;
%                 ks_vaso_new_new =  ks_vaso_new.*exp(1i*(del_phi_xx_v));
                ks_bold_new_new(:,:,i,:) =  ks_bold_new(:,:,i,:).*exp(1i*(del_phi_xx_b));

                
                 %%%%%%%%%%%%%
             end
             ks_vaso_new = ks_vaso_new_new;
             ks_bold_new = ks_bold_new_new;
         end

      %%% Repetition DORK...
        if params.rep_dork == 1
            for i =1:size(ks_vaso,3)
                rep_range = f_nav:l_nav;

                nav_ni_b = (ks_tmp_b(rep_range,(params.gen.n(3)/2)+1,i));
                nav_nn_b = nav.bold(:,(params.gen.n(3)/2)+1,i);
                del_phi_nn_b = (angle(nav_nn_b)-angle(nav_rn_b));%./params.gen.TE; % Nav ph diff
                del_phi_nn_b = mean(del_phi_nn_b);
                del_phi_ni_b = (angle(nav_ni_b)-angle(nav_ri_b));%./params.gen.TE; % Nav ph diff
                del_phi_ni_b = mean(del_phi_ni_b);


                nav_ni_v = (ks_tmp_v(rep_range,(params.gen.n(3)/2)+1,i));
                nav_nn_v = nav.vaso(:,(params.gen.n(3)/2)+1,i);
                del_phi_nn_v = (angle(nav_nn_v)-angle(nav_rn_v));%./params.gen.TE; % Nav ph diff
                del_phi_nn_v = mean(del_phi_nn_v);
                del_phi_ni_v = (angle(nav_ni_v)-angle(nav_ri_v));%./params.gen.TE; % Nav ph diff
                del_phi_ni_v = mean(del_phi_ni_v);


                del_omg_n_b(i) = (del_phi_nn_b-del_phi_ni_b)./(tn-tmp_te);
                del_omg_n_v(i) = (del_phi_nn_v-del_phi_ni_v)./(tn-tmp_te);

                del_phi_n0_b = ((tn*del_phi_ni_b)-(tmp_te*del_phi_nn_b))./(tn-tmp_te);
                del_phi_n0_v = ((tn*del_phi_ni_v)-(tmp_te*del_phi_nn_v))./(tn-tmp_te);

                if params.spi.interl > 1
                    ref_interl = floor(params.spi.interl/2);
                    ref_interl_range = (f_nav)+(params.gen.ro_samples*(ref_interl-1)):(l_nav)+(params.gen.ro_samples*(ref_interl-1));
                    nav0_b_interl = (ks_tmp_b(ref_interl_range,params.gen.n(3)/2+1,i));
                    nav0_v_interl = (ks_tmp_v(ref_interl_range,params.gen.n(3)/2+1,i));
    %                 ks_tmp_b_interl = ks_tmp_b(f_nav:l_nav,:,:);
    %                 ks_tmp_v_interl = ks_tmp_v(f_nav:l_nav,:,:);
                    for j=1:params.spi.interl
                        interl_range = (params.gen.ro_samples*(j-1))+1:params.gen.ro_samples*j;
                        interl_range = interl_range.';
                        nav_range = (f_nav)+(params.gen.ro_samples*(j-1)):(l_nav)+(params.gen.ro_samples*(j-1));
                        nav_b_interl = (ks_tmp_b(nav_range,params.gen.n(3)/2+1,i));
                        del_b_interl = (angle(nav_b_interl)-angle(nav0_b_interl))./tmp_te;
                        del_b_interl = mean(del_b_interl);
                        nav_v_interl = (ks_tmp_v(nav_range,params.gen.n(3)/2+1,i));
                        del_v_interl = (angle(nav_v_interl)-angle(nav0_v_interl))./tmp_te;
                        del_v_interl = mean(del_v_interl);
    %                     % Temp, making it delta negative
    %                     del_b_interl = del_b_interl.*-1;
    %                     del_v_interl = del_v_interl.*-1;
                        ks_bold_new(interl_range,:,i,:) =  ks_bold_new(interl_range,:,i,:).*exp(-1i*del_b_interl.*t(interl_range));
                        ks_vaso_new(interl_range,:,i,:) =  ks_vaso_new(interl_range,:,i,:).*exp(-1i*del_v_interl.*t(interl_range));                 
                    end
                end
    %             del_phi_n0_v = del_phi_n0_v.*-1;
    %             del_phi_n0_b = del_phi_n0_b.*-1;

                ks_bold_new(:,:,i,:) =  ks_bold_new(:,:,i,:).*exp(-1i*(del_phi_n0_b+(del_omg_n_b(i).*t)));
                ks_vaso_new(:,:,i,:) =  ks_vaso_new(:,:,i,:).*exp(-1i*(del_phi_n0_v+(del_omg_n_v(i).*t)));
            end
        end
       
        % Finding Delta W to correctB0 maps
        if params.rep_dork == 1
            %%%%% Trying to get the phase difference to reference after all
            %%%%% corrections, here I assume the reference is 0, since that is
            %%%%% the phase I would expect in the B0 map

            ref_rep = 2;
            ks_tmp_b = ks_bold_new;
            nav_rn_b = squeeze(ks_tmp_b(rep_range,(params.gen.n(3)/2)+1,ref_rep,:));
            nav_rn_b = mean(nav_rn_b,2);
            omg_nav_rn_b = angle(nav_rn_b);%./tmp_te;

            ks_tmp_v = ks_vaso_new;
            nav_rn_v = squeeze(ks_tmp_v(rep_range,(params.gen.n(3)/2)+1,ref_rep,:));
            nav_rn_v = mean(nav_rn_v,2);
            omg_nav_rn_v = angle(nav_rn_v);%./tmp_te;

            for i =1:size(ks_vaso_new,3)
                rep_range = f_nav:l_nav;
                
                nav_ni_b = squeeze(ks_tmp_b(rep_range,(params.gen.n(3)/2)+1,i,:));
                nav_ni_b = mean(nav_ni_b,2);
                omg_nav_ni_b = angle(nav_ni_b);%./tmp_te;

                nav_ni_v = squeeze(ks_tmp_v(rep_range,(params.gen.n(3)/2)+1,i,:));
                nav_ni_v = mean(nav_ni_v,2);
                omg_nav_ni_v = angle(nav_ni_v);%./tmp_te;;


                del_omg_n_b_new(i) = (mean(omg_nav_ni_b-omg_nav_rn_b))./tmp_te;
                del_omg_n_v_new(i) = (mean(omg_nav_ni_v-omg_nav_rn_v))./tmp_te;
                 
%                 del_omg_n_b_new(i) = (mean(omg_nav_ni_b)./tmp_te)-(623.95); % 724.63/ 587
%                 del_omg_n_v_new(i) = (mean(omg_nav_ni_v)./tmp_te)-(623.95);

            end

             save(sprintf('./data/%s/tmp/%s_del_omg_n_v_new.mat',folder,scan),'del_omg_n_v_new')
             save(sprintf('./data/%s/tmp/%s_del_omg_n_b_new.mat',folder,scan),'del_omg_n_b_new')
         
             save(sprintf('./data/%s/tmp/%s_del_omg_n_v.mat',folder,scan),'del_omg_n_v')
             save(sprintf('./data/%s/tmp/%s_del_omg_n_b.mat',folder,scan),'del_omg_n_b')
         
             %%%%%%
        end
       
        % Dirty trick for the moment, to make this part of code work with
        % ABC data...
        if params.gen.seq == 1
            ks_bold = ks_bold_new;
            ks_vaso = ks_vaso_new;
        elseif params.gen.seq == 2
            ks_abc = ks_bold_new;
        end
        clearvars ks_bold_new ks_vaso_new      
end
%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % Permuting to make diemnsions compatible with ISMRM reader in Julia
        ks_traj.kx = permute(ks_traj.kx,[3,2,1]);
        ks_traj.ky = permute(ks_traj.ky,[3,2,1]);
        ks_traj.kz = permute(ks_traj.kz,[3,2,1]);
        
        

        
        if contains(scan,'sv')  ||  contains(scan,'cv')
            ks_vaso = double(ks_vaso);
            ks_bold = double(ks_bold);
        elseif contains(scan,'abc')
            ks_abc = double(ks_abc);
        elseif contains(scan,'sb')
            ks_sb = double(ks_sb);
        end

        % Permuting ks_traj, to make things easier
        kx_tmp = permute(ks_traj.kx,[3 2 1]);
        ky_tmp = permute(ks_traj.ky,[3 2 1]);
        kz_tmp = permute(ks_traj.kz,[3 2 1]);
                
%% Loop over repetitions
% for j = 1:params.repetitions
for j = params.reps_to_save(1):params.reps_to_save(end)
    
                % Checking if dataset exists and delete it
                warning('off')
                if contains(scan,'sv')  ||  contains(scan,'cv')
                    if params.is2d == 1
                        filename_v = sprintf('./data/%s/ismrmd/2d/%s_v_2d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                        filename_b = sprintf('./data/%s/ismrmd/2d/%s_b_2d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                    else
                        filename_v = sprintf('./data/%s/ismrmd/3d/%s_v_3d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                        filename_b = sprintf('./data/%s/ismrmd/3d/%s_b_3d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                    end
                    delete(filename_v);
                    delete(filename_b);
                elseif contains(scan,'abc')
                    if params.is2d == 1
                        filename = sprintf('./data/%s/ismrmd/2d/%s_2d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                    else
                        filename = sprintf('./data/%s/ismrmd/3d/%s_3d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                    end
                    delete(filename);
                elseif contains(scan,'sb')
                    if params.is2d == 1
                        filename = sprintf('./data/%s/ismrmd/2d/%s_2d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                    else
                        filename = sprintf('./data/%s/ismrmd/3d/%s_3d_r%i_%s%s%s.h5',folder,scan,j,traj_name,p_dork,r_dork);
                    end
                    delete(filename);
                end
                warning('on')

                % Create data set and Acquisition
                if contains(scan,'sv')  ||  contains(scan,'cv')
                    dset_v = ismrmrd.Dataset(filename_v);
                    dset_b = ismrmrd.Dataset(filename_b);
                elseif contains(scan,'abc')
                    dset_abc = ismrmrd.Dataset(filename);
                elseif contains(scan,'sb')
                    dset_sb = ismrmrd.Dataset(filename);
                end

                acqData = ismrmrd.Acquisition(params.slices);
                % Set the header elements that don't change
                acqData.head.version(:) = 1;
                acqData.head.number_of_samples(:) = params.nx*params.spi.interl;
                acqData.head.center_sample(:) = floor(params.nx*params.spi.interl/2);
                acqData.head.active_channels(:) = params.ch;
                acqData.head.read_dir  = single(repmat(repmat([1 0 0]',[1 params.ny]),[1 params.slices]));
                acqData.head.phase_dir = single(repmat(repmat([0 1 0]',[1 params.ny]),[1 params.slices]));
                acqData.head.slice_dir = single(repmat(repmat([0 0 1]',[1 params.ny]),[1 params.slices]));
                acqData.head.sample_time_us = single(repmat(params.gen.acqTR,[1 params.slices]));  % In seconds (should be in seconds for MRIReco)

                % AMM:
                if params.is2d == 1
                    acqData.head.trajectory_dimensions(:) = 2;
                else
                    acqData.head.trajectory_dimensions(:) = 3;
                end
            % for k=1:params.slices % slices/partitions
                %% Adding new data
                if contains(scan,'sv')  ||  contains(scan,'cv')

                    % VASO
                    for i=1:last_partition
                      % Set the header elements that change from acquisition to the next
                        % c-style counting
                        acqData.head.scan_counter(i) = params.gen.ro_samples + i-1;
%                         acqData.head.scan_counter(i) = (j-1)*params.gen.ro_samples + i-1;
                        % acqData.head.idx.repetition(i) = j - 1;

                        % Set the flags
                        acqData.head.flagClearAll(i);
                        if i == 1
                            acqData.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_REPETITION', i);
                        elseif i == params.slices
                            acqData.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_LAST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_LAST_IN_REPETITION', i);
                        end  

                        if params.is2d == 1
                            acqData.head.idx.slice(i) = i;
                            traj = [kx_tmp(:,i) ky_tmp(:,i)];
                            data = squeeze(ks_vaso(:,i,j,:));
                        else
                            acqData.head.idx.kspace_encode_step_1(i) = i-1; 
                            traj = [kx_tmp(:,i) ky_tmp(:,i) kz_tmp(:,i)];
                            data = reshape(ks_vaso(:,i,j,:),[],size(ks_vaso,4));
                        end
                        acqData.data{i} = data;
                        acqData.traj{i} = traj.';
                    end
                    
                    dset_v.appendAcquisition(acqData);

                    % BOLD
                    for i=1:last_partition
                        % Set the header elements that change from acquisition to the next
                        % c-style counting
                        acqData.head.scan_counter(i) = params.gen.ro_samples + i-1;
%                         acqData.head.scan_counter(i) = (j-1)*params.gen.ro_samples + i-1;
                        % acqData.head.idx.repetition(i) = j - 1;

                        % Set the flags
                        acqData.head.flagClearAll(i);
                        if i == 1
                            acqData.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_REPETITION', i);
                        elseif i == params.slices
                            acqData.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_LAST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_LAST_IN_REPETITION', i);
                        end  

                        if params.is2d == 1
                            acqData.head.idx.slice(i) = i;
                            traj = [kx_tmp(:,i) ky_tmp(:,i)];
                            data = squeeze(ks_bold(:,i,j,:));
                        else
                            acqData.head.idx.kspace_encode_step_1(i) = i-1; 
                            traj = [kx_tmp(:,i) ky_tmp(:,i) kz_tmp(:,i)];
                            data = reshape(ks_bold(:,i,j,:),[],size(ks_bold,4));
                        end
                        acqData.data{i} = data;
                        acqData.traj{i} = traj.';
                    end
                    
                    dset_b.appendAcquisition(acqData);
                    
                    
                elseif contains(scan,'abc')
                    ks_abc = double(ks_abc);

                    % ABC
                    for i=1:last_partition
                         % Set the header elements that change from acquisition to the next
                        % c-style counting
                        acqData.head.scan_counter(i) = params.gen.ro_samples + i-1;
%                         acqData.head.scan_counter(i) = (j-1)*params.gen.ro_samples + i-1;
                        % acqData.head.idx.repetition(i) = j - 1;

                        % Set the flags
                        acqData.head.flagClearAll(i);
                        if i == 1
                            acqData.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_REPETITION', i);
                        elseif i == params.slices
                            acqData.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_LAST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_LAST_IN_REPETITION', i);
                        end  

                        if params.is2d == 1
                            acqData.head.idx.slice(i) = i;
                            traj = [kx_tmp(:,i) ky_tmp(:,i)];
                            data = squeeze(ks_abc(:,i,j,:));
                        else
                            acqData.head.idx.kspace_encode_step_1(i) = i-1; 
                            traj = [kx_tmp(:,i) ky_tmp(:,i) kz_tmp(:,i)];
                            data = reshape(ks_abc(:,i,j,:),[],size(ks_abc,4));
                        end
                        acqData.data{i} = data;
                        acqData.traj{i} = traj.';
                    end
                    dset_abc.appendAcquisition(acqData);
                    
                elseif contains(scan,'sb')
                    ks_sb = double(ks_sb);

                    % SB
                    for i=1:last_partition
                         % Set the header elements that change from acquisition to the next
                        % c-style counting
                        acqData.head.scan_counter(i) = params.gen.ro_samples + i-1;
%                         acqData.head.scan_counter(i) = (j-1)*params.gen.ro_samples + i-1;
                        % acqData.head.idx.repetition(i) = j - 1;

                        % Set the flags
                        acqData.head.flagClearAll(i);
                        if i == 1
                            acqData.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_FIRST_IN_REPETITION', i);
                        elseif i == params.slices
                            acqData.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i);
                            acqData.head.flagSet('ACQ_LAST_IN_SLICE', i);
                            acqData.head.flagSet('ACQ_LAST_IN_REPETITION', i);
                        end  

                        if params.is2d == 1
                            acqData.head.idx.slice(i) = i;
                            traj = [kx_tmp(:,i) ky_tmp(:,i)];
                            data = squeeze(ks_sb(:,i,j,:));
                        else
                            acqData.head.idx.kspace_encode_step_1(i) = i-1; 
                            traj = [kx_tmp(:,i) ky_tmp(:,i) kz_tmp(:,i)];
                            data = reshape(ks_sb(:,i,j,:),[],size(ks_sb,4));
                        end
                        acqData.data{i} = data;
                        acqData.traj{i} = traj.';
                    end
                    dset_sb.appendAcquisition(acqData);

                end
            % end
        
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
        header.encoding.reconSpace.fieldOfView_mm.x = params.fov(1);
        header.encoding.reconSpace.fieldOfView_mm.y = params.fov(2);
        header.encoding.reconSpace.fieldOfView_mm.z = params.fov(3);
        header.encoding.reconSpace.matrixSize.x = params.mtx_s(1);
        header.encoding.reconSpace.matrixSize.y = params.mtx_s(2);
        header.encoding.reconSpace.matrixSize.z = params.mtx_s(3);
        
        
        % Encoding Limits
        header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
        header.encoding.encodingLimits.kspace_encoding_step_0.maximum = size(traj,1)-1;
        header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(size(traj,1)/2);
        header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
        header.encoding.encodingLimits.kspace_encoding_step_1.maximum = size(traj,2)-1;
        header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(size(traj,3)/2);

        % header.encoding.encodingLimits.repetition.minimum = 0;
        % header.encoding.encodingLimits.repetition.maximum = params.repetitions-1;
        % header.encoding.encodingLimits.repetition.center = 0;

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
        elseif contains(scan,'sb')
            dset_sb.writexml(xmlstring);
            dset_sb.close();  
        end

        fprintf('ISMRMD file, rep %i saved in path %s/ismrmd/ \n',j, folder);
end
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

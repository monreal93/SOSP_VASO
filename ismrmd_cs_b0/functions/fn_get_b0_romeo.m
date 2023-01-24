function fn_get_b0_romeo(folder,scan,file,params)
    
    fprintf('--- Generating B0 maps for scan %s... \n',scan);
    
    cs_mtx = params.mtx_s;
    
    path_save = sprintf('./data/%s/acq/romeo/',folder);
    path_save = what(path_save); path_save = path_save.path;
    path_save = sprintf('%s/b0_%s_%i_%i_%i',path_save,scan,cs_mtx(1),cs_mtx(2),cs_mtx(3));
    save_file_b0 = sprintf('%s/B0_masked_sm.nii',path_save);
    
    tmp = 'y';
    if exist(save_file_b0) > 1
        prompt = 'B0 map file already exists, do you want to replace it: type y/n [y]: ';
        tmp = input(prompt,'s');
    end
    
    if tmp == 'n'
            fprintf('B0 map not generated \n');  
            return
    elseif tmp == 'y'

        var = load(sprintf('./data/%s/raw/%s.mat',folder,file));
        var = struct2cell(var);
        var = var{1};

        sz = size(var);
        x = sz(1);
        y = sz(2);
        z = sz(3);
        ech = sz(4);
        ch = sz(5);
        
%        % Changing resolution of B0 for Sensitivities
%             if cs_mtx(1) > x
%                 var=FFT_1D(var,'kspace',3);
%                 if mod(size(var,3),2)>0
%                     var=padarray(var,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,0],'both');
%                     var=padarray(var,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
%                     var=padarray(var,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
%                 else
%                     var=padarray(var,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,(cs_mtx(3)-z)/2],'both');
%                 end
%                 var=FFT_1D(var,'image',3);
%             else
%                 var = FFT_1D(var,'kspace',3);
%                 if mod(size(var,3),2)>0
%                     var=padarray(var,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
%                     var=padarray(var,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
%                 else
%                     var = padarray(var,[0,0,ceil((cs_mtx(3)-z)/2),0],'both');
%                 end
%                 var = FFT_1D(var,'image',3);
%                 var = FFT_2D(var,'image',1,2);
%                 var = imresize(var,[cs_mtx(1) cs_mtx(2)]);
%                 var = FFT_2D(var,'kspace',1,2);
%         %         var = bart(sprintf('resize -c 0 %i 1 %i',cs_mtx(1),cs_mtx(2)),var);
%             end
            

        img = FFT_2D(var,'image',1,2);
        % Shifting slices, not sure if always needed
        % img = circshift(img,-1,3);

        % Creating romeo Directory
        path_read = sprintf('./data/%s/acq/romeo',folder);
        path_read = what(path_read); path_read = path_read.path;
        
        warning('off')
        % system(sprintf('mkdir %s',path_read));
        % Creating specific folder for mtx size
        system(sprintf('mkdir %s',path_save));
        warning('on')
        
        % Subsetting if is 2D
        if cs_mtx(3) == 1
            % Writing Niftis of mag and phase
            img_ph = angle(img);
            img_mag = abs(img);
            niftiwrite(img_mag(:,:,z/2+1,:),sprintf('%s/%s_mag.nii',path_save,'fieldmap'))
            niftiwrite(img_ph(:,:,z/2+1,:),sprintf('%s/%s_ph.nii',path_save,'fieldmap'))
        else        
            % Writing Niftis of mag and phase
            img_ph = angle(img);
            img_mag = abs(img);
            niftiwrite(img_mag,sprintf('%s/%s_mag.nii',path_save,'fieldmap'))
            niftiwrite(img_ph,sprintf('%s/%s_ph.nii',path_save,'fieldmap'))
        end
     
        % Writing resized Ech1 and Ech9 for future refererence
        tmp = rssq(img_mag,5);
        % tmp = imresize(tmp,[params.gen.n(1) params.gen.n(2)]);
       % Changing resolution of echos to save reference
        if cs_mtx(1) > x
            tmp=FFT_3D(tmp,'kspace');
            if mod(size(tmp,3),2)>0
                tmp=padarray(tmp,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,0],'both');
                tmp=padarray(tmp,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
                tmp=padarray(tmp,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
            else
                if cs_mtx(3) == 1
                    tmp=padarray(tmp,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,0],'both');
                else
                    tmp=padarray(tmp,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,(cs_mtx(3)-z)/2],'both');
                end
            end
            tmp=FFT_3D(tmp,'image');
        else
            tmp = FFT_1D(tmp,'kspace',3);
            if mod(size(tmp,3),2)>0
                tmp=padarray(tmp,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
                tmp=padarray(tmp,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
            else
                tmp = padarray(tmp,[0,0,ceil((cs_mtx(3)-z)/2),0],'both');
            end
            tmp = FFT_1D(tmp,'image',3);
        end
        % Subsetting if is 2D
        if cs_mtx(3) == 1
                    tmp = tmp(:,:,z/2+1,:);
        end
            
        % I am flipping dim 1 just to match recon:
        ech1 = flip(tmp(:,:,:,1),1)*1e5;
        ech9 = flip(tmp(:,:,:,9),1)*1e5;
        niftiwrite(abs(ech1),sprintf('%s/gre_%i_%i_%i_ech1.nii',path_save,params.gen.n(1),params.gen.n(2),params.gen.n(3)))
        niftiwrite(abs(ech9),sprintf('%s/gre_%i_%i_%i_ech9.nii',path_save,params.gen.n(1),params.gen.n(2),params.gen.n(3)))

        % Calling ROMEO
        % tmp = '[8.14,17.11,23,25.01,27.02,32.11,55,60,65]';  % This has to be read from the .dat file and saved in a variable
        % For 9 echos fieldmap
        te = string((params.twix_params_b0.TE(1:9))/1000);
       te = strcat('[',te(1),',',te(2),',',te(3),',',te(4),',',te(5),',',te(6),',',te(7),',',te(8),',',te(9),']');
        % For 3 echos fieldmap
%         te = string((params.twix_params_b0.TE(1:4))/1000);
%         te = strcat('[',te(1),',',te(2),',',te(3),',',te(4),']');
        romeo_settings = ...
            sprintf('romeo -p %s/%s_ph.nii -m %s/%s_mag.nii -k nomask -t %s -o %s --compute-B0 --phase-offset-correction', ...
            path_save,file,path_save,file,te,path_save);
        
        % te = (params.twix_params_b0.TE(1:9))/1000;
        parameters.output_dir = path_save;
        parameters.TE = te;
        parameters.mag = load_nii(sprintf('%s/%s_mag.nii',path_save,'fieldmap'));
        parameters.mag = parameters.mag.img;
        parameters.mask = 'nomask';
        parameters.calculate_B0 = true;
        parameters.phase_offset_correction = 'off';
        parameters.voxel_size = load_nii_hdr(sprintf('%s/%s_ph.nii',path_save,'fieldmap'));
        parameters.voxel_size = parameters.voxel_size.dime.pixdim(2:4);
        parameters.additional_flags = '--verbose -q -i';
        phase = load_nii(sprintf('%s/%s_ph.nii',path_save,'fieldmap'));
        phase = phase.img;
        
        fprintf('Calculating the B0 maps using ROMEO...')
        % [unwrapped,b0] = ROMEO(phase,parameters);
        % system(romeo_settings);
        fprintf('B0 map saved as path %s/B0.nii \n',path_save);
        
        fprintf('For now ROMEO from Matlab is not working.... paste this in the command line and run script again')
        fprintf(sprintf('\n %s', romeo_settings)) 

        % Masking
        b0 = niftiread(sprintf('%s/B0.nii',path_save));
        msk = rssq(img_mag,5); msk = msk(:,:,:,1);
        msk(msk>1e-4) = 1;
        msk(msk<1e-4) = 0;
        b0 = msk.*real(b0);
        b0(isnan(b0)) = 0;
        
       % Changing resolution of B0 for Sensitivities
%         if cs_mtx(1) > x
%             b0=FFT_3D(b0,'kspace');
%             if mod(size(b0,3),2)>0
%                 b0=padarray(b0,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,0],'both');
%                 b0=padarray(b0,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
%                 b0=padarray(b0,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
%             else
%                 b0=padarray(b0,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,(cs_mtx(3)-z)/2],'both');
%             end
%             b0=FFT_3D(b0,'image');
%         else
%             b0 = FFT_1D(b0,'kspace',3);
%             if mod(size(b0,3),2)>0
%                 b0=padarray(b0,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
%                 b0=padarray(b0,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
%             else
%                 b0 = padarray(b0,[0,0,ceil((cs_mtx(3)-z)/2),0],'both');
%             end
%             b0 = FFT_1D(b0,'image',3);
% %             b0 = FFT_2D(b0,'image',1,2);
% %             b0 = imresize(b0,[cs_mtx(1) cs_mtx(2)]);
%         end
        b0 = imresize(b0,[cs_mtx(1) cs_mtx(2)]);
        b0 = real(b0);
        
        % Saving masked map
        niftiwrite(b0,sprintf('%s/B0_masked.nii',path_save));

        % Smoothing, for now just a gauss filter
        % b0 = imgaussfilt(b0,2);
        % Remove gaps and smooth
        tmp = abs(ech1);
        for s=1:z
            disp(['Sli #' num2str(s)])
            b0(:,:,s) = FitB0ToRemoveGaps(b0(:,:,s),tmp(:,:,s),5);
        end

        % Saving masked smoothed map
        niftiwrite(b0,sprintf('%s/B0_sm.nii',path_save));

        % Saving mask
        niftiwrite(msk,sprintf('%s/mask.nii',path_save));
    end
end

% Script to read and format twix data



%% Select the name of file...
% psfy_ref, psfy_on, psfz_ref, psfz_on
% wc, non_wc, u_wc, u_non_wc
% fieldmap
% sv_1p2, sv_0p8, epi_1p2, epi_0p9 , b0

function [twix_params,twix_params_b0] = fn_read_twix(folder,scan,params)

    fprintf('--- Reading scan %s Twix data and saving it into .mat format... \n',scan);
    tmp = 'y';
    twix_params_b0 = [];
    twix_params = [];
    
    % Here I need to list the different types of seq I will use
    if  contains(scan,'sv') ||  contains(scan,'cv')
        save_file_vaso = sprintf('./data/%s/raw/%s_ks_vaso.mat',folder,scan);
        save_file_bold = sprintf('./data/%s/raw/%s_ks_bold.mat',folder,scan);
    elseif contains(scan,'b0')
        save_file = sprintf('./data/%s/raw/%s.mat',folder,scan);
    elseif contains(scan,'abc')
        save_file = sprintf('./data/%s/raw/%s_ks_abc.mat',folder,scan);
    elseif contains(scan,'sb')
        save_file = sprintf('./data/%s/raw/%s_ks_sb.mat',folder,scan);
    end
    
    if contains(scan,'sv') || contains(scan,'cv')
        if exist(save_file_vaso,'file') > 1
            prompt = 'File already exists, do you want to replace it: type y/n [y]: ';
            tmp = input(prompt,'s');
        end
    elseif contains(scan,'b0') || contains(scan,'abc') || contains(scan,'sb')
        if exist(save_file,'file') > 1
            prompt = 'File already exists, do you want to replace it: type y/n [y]: ';
            tmp = input(prompt,'s');
        end
    end
        if tmp == 'n'
            fprintf('Formated .mat twix data not generated \n');  
            return
        elseif tmp == 'y'
            r = params.rz;
            pf = params.pf;
%             sl = floor(params.slices/r*pf);
            sl = params.mtx_s(3);
            
            ro = params.nx;
            ch = params.ch;
            interl = params.spi.interl;

            % twix = mapVBVD();
%             twix = mapVBVD(sprintf('./data/%s/raw/twix/%s.dat',folder,scan));
            % Find the file that matches the name of scan
            tmp = dir(sprintf('./data/%s/raw/twix/*%s*',folder,scan));
            twix = mapVBVD(sprintf('%s/%s',tmp.folder,tmp.name));
%             twix = twix{1}; % AMM: This seems to be only needed for VE12...

            % Saving relevant parameters
            if contains(scan,'b0') 
                % Here I can save more twix parameters as needed...
                if contains(folder,'2023') && params.gen.field_strength == 7
                    twix = twix{1};
                end
               
                twix_params_b0.shift(1) = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor; 
                twix_params_b0.shift(2) = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra; 
                twix_params_b0.TE = twix.hdr.Meas.alTE(1:5);
                twix_params_b0.TR = twix.hdr.Meas.alTR;
                twix_params_b0.FA = twix.hdr.Phoenix.adFlipAngleDegree;
                twix_params_b0.scan_time = twix.hdr.Phoenix.lScanTimeSec;
                twix_params_b0.ph_fov = twix.hdr.Meas.PeFOV;
                twix_params_b0.dwell = twix.hdr.Meas.alDwellTime;
                % Calculating the scan resolution
                twix_params_b0.res = [twix.hdr.Meas.RoFOV/twix.hdr.Meas.BaseResolution ...
                    twix.hdr.Meas.PhaseFoV/twix.hdr.Meas.PhaseEncodingLines ...
                    twix.hdr.Meas.dSliceResolution].*1e-3;
            else
                % Here I can save more twix parameters as needed...
                % the parameter where version is saved is
                % twix.image.softwareVersion... 
                % For VE12
                if contains(folder,'2023') && params.gen.field_strength == 7
                    twix = twix{1};
                end
                
                % First I make the shift zero, then I update, need to check
                % why in 9.4 I dont se it or tilted?
                params.gen.shift(1:2) = 0;
                twix_params.shift(1) = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor; 
                twix_params.shift(2) = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra;
                twix_params.TE = twix.hdr.Meas.alTE(1).*1e-6;
                twix_params.TR = twix.hdr.Meas.alTR(1).*1e-6;
                twix_params.FA = twix.hdr.Phoenix.adFlipAngleDegree{1};
                twix_params.scan_time = twix.hdr.Phoenix.lScanTimeSec;
                twix_params.ph_fov = twix.hdr.Meas.PeFOV.*1e-3;
                twix_params.dwell = twix.hdr.Meas.alDwellTime(1).*1e-10; % not sure about this units
                twix_params.mtx_s = params.mtx_s;
                twix_params.ch = params.ch;
                twix_params.repetitions = params.repetitions;
%                 twix_params_sv.sl_to_recon = params.slice_to_save;
                twix_params.rz = params.rz;
                twix_params.grad_delay = [twix.hdr.Dicom.lGradDelayTimeX, ...
                    twix.hdr.Dicom.lGradDelayTimeY, twix.hdr.Dicom.lGradDelayTimeZ]*1e-6;
            end

            if contains(scan,'b0')
                twix.image.flagRemoveOS = 1;
            else
                twix.image.flagRemoveOS = 0;
            end
            twix.image.flagDoAverage = 1;

            % Get slice info
            read = twix.hdr.Config.NImageCols;
            p_enc = twix.hdr.Config.NImageLins/r;
            avg = twix.hdr.Config.NAve;
            % ch = twix.hdr.Config.NChaMeas;

            if contains(scan,'b0')
                image_data = twix.image(); 
                 aa=squeeze(image_data);
                %  aa=mean(aa,5);  % Averaging echos
            %      aa=aa(:,:,:,:,3); % Taking one echo
                 aa=permute(aa,[1,3,4,5,2]);      % Moving Ch dim to end
                 % Re-ordering slices due to interleaved acquisition
                 bb=zeros(size(aa));
                 bb(:,:,2:2:end,:)=aa(:,:,1:floor(size(aa,3)/2),:);
                 bb(:,:,1:2:end,:)=aa(:,:,floor(size(aa,3)/2)+1:end,:);
            %      bb=FFT_1D(bb,'kspace',3);
            %      bb=padarray(bb,[52,64,6],'both');
            %      bb=FFT_1D(bb,'image',3);
                 b0=bb;
                 % AMM: Not sure if this is right.. Double check...
                % Here I do FFT Shift to try to make it in the isocenter
%                 shx = (twix_params_b0.shift*1e-4)./twix_params_b0.res(1);
%                 b0 = fft_shift_2D(b0,1,shx,1,0);
            elseif contains(scan,'sv')
                dyn = (params.repetitions)*2;       % Here I do *2 because I have VASO-BOLD
                image_data = twix.image();
                % trying to reshare matrix, for Spiral acquisitions
                aa=squeeze(image_data);
                % AMM: This only applies if I use sync scan, now I don't
                % have it so I dont need to remove skope sync scans.
%                 if params.gen.skope
%                     bb=aa(:,:,sl+1:end,:);  % Removing the skope sync readouts
%                 else
                    bb=aa;
%                 end


                cc=permute(bb,[1 4 3 2]); % Getting ch dim to the end and set to second
                
%                 Temp, trying to extract the dork navigators:
                if params.gen.dork == 2
                 if size(cc,3) > sl*dyn*interl
                    dork_nav = cc(:,:,2:2:end,:);
                    dork_nav = reshape(dork_nav,[],sl*dyn,ch);
                    dork_nav = dork_nav(params.gen.ro_samples-99:params.gen.ro_samples,:,:);
                    dork_nav = squeeze(dork_nav);
                    dork_nav = permute(dork_nav,[2 1 3]);
                    dork_nav = reshape(dork_nav,sl*2,params.repetitions,[],ch);
                    dork_nav = permute(dork_nav,[3 1 2 4]);
                    nav.vaso = dork_nav(:,1:sl,:,:);
                    nav.bold = dork_nav(:,sl+1:end,:,:);
                    cc = cc(:,:,1:2:end,:);
                    save(sprintf('./data/%s/raw/%s_nav.mat',folder,scan),'nav')
                 end
                end
                
                dd=reshape(cc,[],sl*dyn,ch);  % Permuting first and second dim, Col and Set

                % VASO contrast
                j = 1;
                for i=1:sl:(dyn/2)*sl %Number of dynamics/2
                    ks_vaso(:,i:i+(sl-1),:) = dd(:,j:j+(sl-1),:);
                    j = j+(sl*2);
                end

                % BOLD contrast
                 j = sl+1;
                for i=1:sl:(dyn/2)*sl %Number of dynamics/2
                    ks_bold(:,i:i+(sl-1),:) = dd(:,j:j+(sl-1),:);
                    j = j+(sl*2);
                end
		
                    % Extracting the dynamics
                    ks_vaso = permute(ks_vaso,[1,3,2]);
                    ks_vaso = reshape(ks_vaso,[ro*interl ch sl dyn/2]);
                    ks_vaso = permute(ks_vaso,[1,3,4,2]);

                    ks_bold = permute(ks_bold,[1,3,2]);
                    ks_bold = reshape(ks_bold,[ro*interl ch sl dyn/2]);
                    ks_bold = permute(ks_bold,[1,3,4,2]);
            elseif contains(scan,'abc')
                dyn = params.repetitions;       % Here I do *2 because I have VASO-BOLD
                image_data = twix.image();
                % trying to reshare matrix, for Spiral acquisitions
                aa=squeeze(image_data);
                % AMM: This only applies if I use sync scan, now I don't
                % have it so I dont need to remove skope sync scans.
%                 if params.gen.skope
%                     bb=aa(:,:,sl+1:end,:);  % Removing the skope sync readouts
%                 else
                    bb=aa;
%                 end
                cc=permute(bb,[1 4 3 2]); % Getting ch dim to the end and set to second
                dd=reshape(cc,[],sl*dyn,ch);  % Permuting first and second dim, Col and Set
		
                % Extracting the dynamics
                ks_abc = permute(dd,[1,3,2]);
                ks_abc = reshape(ks_abc,[ro*interl ch sl dyn]);
                ks_abc = permute(ks_abc,[1,3,4,2]);
                
            elseif contains(scan,'sb')
                dyn = params.repetitions;       % Here I do *2 because I have VASO-BOLD
                image_data = twix.image();
                % trying to reshare matrix, for Spiral acquisitions
                aa=squeeze(image_data);
                % AMM: This only applies if I use sync scan, now I don't
                % have it so I dont need to remove skope sync scans.
%                 if params.gen.skope
%                     bb=aa(:,:,sl+1:end,:);  % Removing the skope sync readouts
%                 else
                    bb=aa;
%                 end
                cc=permute(bb,[1 4 3 2]); % Getting ch dim to the end and set to second
                dd=reshape(cc,[],sl*dyn,ch);  % Permuting first and second dim, Col and Set
		
                % Extracting the dynamics
                ks_sb = permute(dd,[1,3,2]);
                ks_sb = reshape(ks_sb,[ro*interl ch sl dyn]);
                ks_sb = permute(ks_sb,[1,3,4,2]);
                
            elseif contains(scan,'cv')
                dyn = (params.repetitions)*2;       % Here I do *2 because I have VASO-BOLD
                image_data = twix.image();
                % trying to reshare matrix, for Spiral acquisitions
                aa=squeeze(image_data);
                % AMM: This only applies if I use sync scan, now I don't
                % have it so I dont need to remove skope sync scans.
%                 if params.gen.skope
%                     bb=aa(:,:,sl+1:end,:);  % Removing the skope sync readouts
%                 else
                    bb=permute(aa,[3 1 2]);
%                 end
                
%                 und_dim = round(size(bb,1)/params.repetitions/2/params.slices);
                und_dim = size(bb,1)/params.repetitions/params.slices/2;
                und_dim = und_dim-3;
                ks = zeros(size(bb,2),und_dim,params.slices*2,params.repetitions,params.ch);
                nav = zeros(size(bb,2),3,params.slices*2,params.repetitions,params.ch);
                
                tmp = 0;
                for i=1:params.repetitions
                    for j=1:params.slices*2
                        % if mod(j,2) == 1
                            tmp1 = tmp;
                            for k=1:3
                                tmp1 = tmp1+1;
                                nav(:,k,j,i,:) = permute(bb(tmp1,:,:),[2 1 4 5 3]);             
                            end
                            tmp = tmp+3;   
                        % end
                        for k=1:und_dim
                            tmp = tmp + 1;
                            ks(:,k,j,i,:) = permute(bb(tmp,:,:),[2 1 4 5 3]);             
                        end
%                         tmp = j*und_dim-1;
                    end
                    % tmp = i*j*und_dim;
                end
                
                % I am only getting the navigators of Vaso acq
                % nav = nav(:,:,1:2:end,:,:);
                ks_new = ks;
                
%                 % Zero padding...
%                 ks_new = zeros(size(ks));
%                 ks_new = padarray(ks_new,[0 (und_dim*params.epi.ry-und_dim) 0 0 0],'pre');
%                 ks_new(:,1:params.epi.ry:end,:,:,:) = ks;
%                 ks_new = padarray(ks_new,[0 params.gen.n(2)-size(ks_new,2) 0 0 0],'pre'); 
                
%                 cc=permute(bb,[3 1 2]); % Getting ch dim to the end and set to second
%                 
%                 dd=reshape(cc,200,sl*dyn,ro,ch);  % Permuting first and second dim, Col and Set
% %                 dd = permute(dd,[1 3 2 4]);
%                 
                ks_vaso = ks_new(:,:,1:sl,:,:);
                ks_bold = ks_new(:,:,sl+1:end,:,:);
                
            end


            if  contains(scan,'sv')
%                 fn_create_ismrmd5(folder,scan,params,ks_vaso,ks_bold);
                save(save_file_vaso,'ks_vaso','-v7.3')
                save(save_file_bold,'ks_bold','-v7.3')
                save(sprintf('./data/%s/acq/%s_twix_params.mat',folder,scan),'twix_params');
                fprintf('Spiral/Cartesian Vaso data saved in path %s/raw/ \n',folder);
            elseif  contains(scan,'cv')
                save(sprintf('./data/%s/raw/%s_nav.mat',folder,scan),'nav')
                save(save_file_vaso,'ks_vaso','-v7.3')
                save(save_file_bold,'ks_bold','-v7.3')
                save(sprintf('./data/%s/acq/%s_twix_params.mat',folder,scan),'twix_params');
                fprintf('Cartesian Vaso data saved in path %s/raw/ \n',folder);
            elseif contains(scan,'b0')
                save(save_file,'b0','-v7.3')
                save(sprintf('./data/%s/acq/twix_params_b0.mat',folder),'twix_params_b0');
                fprintf('B0 scan saved in path %s/raw/ \n',folder);
            elseif contains(scan,'abc')
                save(save_file,'ks_abc','-v7.3')
                save(sprintf('./data/%s/acq/%s_twix_params.mat',folder,scan),'twix_params');
                fprintf('ABC scan saved in path %s/raw/ \n',folder);
            elseif contains(scan,'sb')
                save(save_file,'ks_sb','-v7.3')
                save(sprintf('./data/%s/acq/%s_twix_params.mat',folder,scan),'twix_params');
                fprintf('Spiral BOLD scan saved in path %s/raw/ \n',folder);
            end
        end 
end

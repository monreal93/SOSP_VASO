% Script to read and format twix data



%% Select the name of file...
% psfy_ref, psfy_on, psfz_ref, psfz_on
% wc, non_wc, u_wc, u_non_wc
% fieldmap
% sv_1p2, sv_0p8, epi_1p2, epi_0p9 , b0

function [twix_params_sv,twix_params_b0] = fn_read_twix(folder,scan,params)

    fprintf('--- Reading scan %s Twix data and saving it into .mat format... \n',scan);
    tmp = 'y';
    twix_params_b0 = [];
    twix_params_sv = [];
    if  contains(scan,'sv')
        save_file_vaso = sprintf('./data/%s/raw/%s_ks_vaso.mat',folder,scan);
        save_file_bold = sprintf('./data/%s/raw/%s_ks_bold.mat',folder,scan);
    elseif contains(scan,'b0')
        save_file = sprintf('./data/%s/raw/%s.mat',folder,scan);
    end
    
    if contains(scan,'sv')
        if exist(save_file_vaso,'file') > 1
            prompt = 'File already exists, do you want to replace it: type y/n [y]: ';
            tmp = input(prompt,'s');
        end
    elseif contains(scan,'b0')
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
            sl = params.slices;
            dyn = (params.repetitions)*2;       % Here I do *2 because I have VASO-BOLD
            ro = params.nx;
            ch = params.ch;
            interl = params.spi.interl;

            % twix = mapVBVD();
%             twix = mapVBVD(sprintf('./data/%s/raw/twix/%s.dat',folder,scan));
            % Find the file that matches the name of scan
            tmp = dir(sprintf('./data/%s/raw/twix/*%s*',folder,scan));
            twix = mapVBVD(sprintf('%s/%s',tmp.folder,tmp.name));

            % Saving relevant parameters
            if contains(scan,'b0')
                % Here I can save more twix parameters as needed...
                twix_params_b0.shift = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor; 
%                 twix_params_b0.shift = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra; 
                twix_params_b0.TE = twix.hdr.Meas.alTE;
                twix_params_b0.TR = twix.hdr.Meas.alTR;
                twix_params_b0.FA = twix.hdr.Phoenix.adFlipAngleDegree;
                twix_params_b0.scan_time = twix.hdr.Phoenix.lScanTimeSec;
                twix_params_b0.ph_fov = twix.hdr.Meas.PeFOV;
                twix_params_b0.dwell = twix.hdr.Meas.alDwellTime;
            else
                % Here I can save more twix parameters as needed...
                twix_params_sv.TE = twix.hdr.Meas.alTE(1).*1e-6;
                twix_params_sv.TR = twix.hdr.Meas.alTR(1).*1e-6;
                twix_params_sv.FA = twix.hdr.Phoenix.adFlipAngleDegree{1};
                twix_params_sv.scan_time = twix.hdr.Phoenix.lScanTimeSec;
                twix_params_sv.ph_fov = twix.hdr.Meas.PeFOV.*1e-3;
                twix_params_sv.dwell = twix.hdr.Meas.alDwellTime(1).*1e-10; % not sure about this units
                twix_params_sv.mtx_s = params.mtx_s;
                twix_params_sv.ch = params.ch;
                twix_params_sv.repetitions = params.repetitions;
%                 twix_params_sv.sl_to_recon = params.slice_to_save;
                twix_params_sv.rz = params.rz;
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

            if contains(scan,'fieldmap')

                sli = twix.hdr.Config.NImagePar/r;
                image_data = twix.image();
                img = squeeze(image_data);
                img = img(:,:,:,:,1);

                img = permute(img,[1 3 4 2]);

                img = padarray(img,[88 88 6 0],'both');

                fieldmap = FFT_3D(img,'image');
            elseif contains(scan,'psf')
                img = squeeze(twix.image());
                    if avg > 1 && ch ==1
                        %img = img;
                    else
                        img = permute(img,[1 3 2]);
                    end
            elseif contains(scan,'b0')

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
                %  aa=FFT_2D(aa,'image',1,2);
            elseif contains(scan,'sv')
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
            end

            if  contains(scan,'sv')
                    save(save_file_vaso,'ks_vaso','-v7.3')
                    save(save_file_bold,'ks_bold','-v7.3')
                    save(sprintf('./data/%s/acq/%s_twix_params.mat',folder,scan),'twix_params_sv');
                    fprintf('Spiral Vaso data saved in path %s/raw/ \n',folder);
            elseif contains(scan,'b0')
                save(save_file,'b0','-v7.3')
                save(sprintf('./data/%s/acq/twix_params_b0.mat',folder),'twix_params_b0');
                fprintf('B0 scan saved in path %s/raw/ \n',folder);
            end
        end 
end

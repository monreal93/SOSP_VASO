function fn_coil_sensitivities(folder,scan,file,params)

    fprintf('--- Generating Sensitivity maps for scan %s... \n',scan);

    cs_mtx = params.mtx_s;
    ncc = params.ncc;
    tmp = 'y';
    
    if ncc > 0
        save_file = sprintf('./data/%s/acq/cs_%i_%i_%i_cc%i.mat',folder,cs_mtx(1),cs_mtx(2),cs_mtx(3),ncc);
    else
        save_file = sprintf('./data/%s/acq/cs_%i_%i_%i.mat',folder,cs_mtx(1),cs_mtx(2),cs_mtx(3));
    end
    
    if exist(save_file,'file') > 1
        prompt = 'Sensitivy Maps file already exists, do you want to replace it: type y/n [y]: ';
        tmp = input(prompt,'s');
    end   
        if tmp == 'n'
            fprintf('Coil Sensitivities not generated \n');  
            return
        elseif tmp == 'y'
          var = load(sprintf('./data/%s/raw/%s.mat',folder,file));
            var = struct2cell(var);
            var = var{1};

            % % Temp
            % var = var(:,:,(24/2)-3:(24/2)+2,:);

            sz = size(var);
            x = sz(1);
            y = sz(2);
            z = sz(3);
            ech = sz(4);
            ch = sz(5);
            slwin = 5; % sliding window length for GCC

            % Taking only one echo
            var = squeeze(var(:,:,:,1,:));

            % Changing resolution of B0 for Sensitivities
            if cs_mtx(1) > x
                var=FFT_1D(var,'kspace',3);
                if mod(size(var,3),2)>0
                    var=padarray(var,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,0],'both');
                    var=padarray(var,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
                    var=padarray(var,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
                else
                    var=padarray(var,[(cs_mtx(1)-x)/2,(cs_mtx(2)-y)/2,(cs_mtx(3)-z)/2],'both');
                end
                var=FFT_1D(var,'image',3);
            else
                var = FFT_1D(var,'kspace',3);
                if mod(size(var,3),2)>0
                    var=padarray(var,[0,0,floor((cs_mtx(3)-z)/2)],'pre');
                    var=padarray(var,[0,0,(floor((cs_mtx(3)-z)/2))+1],'post');
                else
                    var = padarray(var,[0,0,ceil((cs_mtx(3)-z)/2),0],'both');
                end
                var = FFT_1D(var,'image',3);
                var = FFT_2D(var,'image',1,2);
                var = imresize(var,[cs_mtx(1) cs_mtx(2)]);
                var = FFT_2D(var,'kspace',1,2);
        %         var = bart(sprintf('resize -c 0 %i 1 %i',cs_mtx(1),cs_mtx(2)),var);
            end

            % non_wc_crop = non_wc((x/2)-(y/2):(x/2)+(y/2)-1,:,:,:);
            msk = zeros(size(var));
            id = rssq(var,4)>1e-4;
            id = repmat(id,1,1,1,ech,1);
            msk(id) = 1;

            % ks = FFT_2D(var,'kspace',1,2);
            ks=var;
            sensemaps = zeros(size(ks));

            if ncc > 0
                espiritmaps = zeros(size(var,1),size(var,2),size(var,3),ncc);
            else
                espiritmaps = zeros([cs_mtx ch]);
            end

            %% Coil Compression
            if ncc > 0
            %     coil_s = FFT_3D(CoilSensitivity,'kspace');
                ncalib = 24; % use 24 calibration lines to compute compression
                dim = 2;
                [sx,sy,sz,Nc] = size(ks);
                calib = crop(ks,[ncalib,sy,sz,Nc]);
            %     eccmtx = calcECCMtx(calib,dim,ncc);
                gccmtx = calcGCCMtx(calib,dim,slwin);
                % crop and align matrices 
            %     eccmtx_aligned = alignCCMtx(eccmtx(:,1:ncc,:));
                gccmtx_aligned = alignCCMtx(gccmtx(:,1:ncc,:));
                % compress the data
            %     ECCDATA_aligned = CC(ks,eccmtx_aligned, dim);
                CCDATA_aligned = CC(ks,gccmtx_aligned, dim);
            %     ks = ECCDATA_aligned;
                ks = CCDATA_aligned;
                non_wc = FFT_2D(ks,'image',1,2);
            end

            % msk1 = zeros(size(msk));
            % msk1((y/2)-10:(y/2)+9,(y/2)-10:(y/2)+9,:,:) = 1;
            % 
            % ks = ks.*msk1;

            % Generating coil Sensitivities
            for i=1:size(var,3)
            %     sensemaps(:,:,i,:) = bart('caldir 20', ks(:,:,i,:));
            %     espiritmaps(:,:,i,:,:) = bart('ecalib -r 20 -m 2', ks(:,:,i,:));
                espiritmaps(:,:,i,:) = bart('ecalib -r 60 -m 1', ks(:,:,i,:));
            end

            coil_sens = double(squeeze(espiritmaps(:,:,:,:,1)));
            % sens = sensemaps.*msk;
            % espi = espiritmaps.*msk;
            if ncc > 0
                save(save_file,'coil_sens','-v7.3')
            else
                save(save_file,'coil_sens','-v7.3')
            end

            fprintf('Coil Sensitivities saved in path %s/acq/ \n',folder);  
        end

end

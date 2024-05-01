clear all; clc

cd /mnt/5T4/Alejandro/sosp_vaso/
addpath(genpath('../tools/as'))

folder = '08232023_sv_abc_paper';
scan = 'cv_51';
b0 = '';           % b0='','fessler','romeo','gilad','skope','fessler'
traj = 'nom';           % 'nom', 'sk'
rep_recon = 1:40;          % Range or repetitions to merge
is2d = 0;
contrasts = ["b","v"];  % "v","b","abc"
p_dork = "";              % Partition DORK "_pDORK"
r_dork = "_rDORK";              % Repetition DORK "_rDORK"
drift = "";               % Drift "_drift"
rotate = false;

% root = '/mnt/ssh/var/www/dabeast/5T3/Alejandro/sosp_vaso';
root = '/mnt/5T4/Alejandro/sosp_vaso';

if contains(scan,'sv') || contains(scan,'abc') || contains(scan,'sb')
    load(sprintf('%s/data/%s/acq/%s_params.mat',root,folder,scan))

    if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
    %     contrasts = ["v","b"];
        merged = zeros([params.gen.n length(rep_recon)*2]);
    else
    %     contrasts = "abc";
        merged = zeros([params.gen.n length(rep_recon)]);
    end

    k = 1;
    for j=rep_recon(1):rep_recon(end)
        if is2d == 1
            for i = 1:params.gen.n(3)
                    if  isempty(b0)
                        tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s%s%s%s.nii',root,folder,scan,contrasts(1),i,j,traj,p_dork,r_dork,drift));
                        merged(:,:,i,k) = tmp;
                        if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                            tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s%s%s%s.nii',root,folder,scan,contrasts(2),i,j,traj,p_dork,r_dork,drift));
                            merged(:,:,i,k+1) = tmp;
                        end
                    else
                        tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s%s%s%s_b0_%s.nii',root,folder,scan,contrasts(1),i,j,traj,p_dork,r_dork,drift,b0));
                        merged(:,:,i,k) = tmp;
                        if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                            tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s%s%s%s_b0_%s.nii',root,folder,scan,contrasts(2),i,j,traj,p_dork,r_dork, drift,b0));
                            merged(:,:,i,k+1) = tmp;
                        end
                    end
            end
        else
            if  isempty(b0)
                    tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s%s%s%s.nii',root,folder,scan,contrasts(1),j,traj,p_dork,r_dork,drift));
                    merged(:,:,:,k) = tmp;
                    if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                        tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s%s%s%s.nii',root,folder,scan,contrasts(2),j,traj,p_dork,r_dork,drift));
                        merged(:,:,:,k+1) = tmp;
                    end
            else
                    tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s%s%s%s_b0_%s.nii',root,folder,scan,contrasts(1),j,traj,p_dork,r_dork,drift,b0));
    %                 tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s_b0_%s_mrreco_noDORK.nii',root,folder,scan,contrasts(1),j,traj,b0));
                    merged(:,:,:,k) = tmp;
                    if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                        tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s%s%s%s_b0_%s.nii',root,folder,scan,contrasts(2),j,traj,p_dork,r_dork,drift,b0));
                        merged(:,:,:,k+1) = tmp;
                    end
            end
        end
    if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
        k = k+2;
    else
        k = k+1;
    end
    end

    merged = merged.*1e4;

    % % Sometimes I need to scale it even more
    merged = merged.*1e6;

    % Adding nifti information
    % info = niftiinfo(sprintf('./data/%s/recon/3d/%s_v_rep%i_3d_mrreco.nii',folder,scan,1));
    % Reading dummy nifti to get info
    info =  niftiinfo('brain.nii');

    % AMM: Need to double check this parameters:
    info.Description = sprintf('Scan %s %s',folder,scan);
    info.SliceCode = 'Sequential-Increasing';
    %%%%
    info.PixelDimensions = params.gen.res*1000;
    info.ImageSize = params.gen.n;
    info.ImageSize(3) = info.ImageSize(3).*params.gen.kz;
    if length(rep_recon) > 1
        info.PixelDimensions(4) = params.gen.volTR;
        if length(rep_recon)>1 && length(contrasts) >1
            info.ImageSize(4) = length(rep_recon)*2;
        else
            info.ImageSize(4) = length(rep_recon);
        end
    end
    info.SpaceUnits = 'Millimeter';
    info.TimeUnits = 'Second';
    info.Datatype = 'single';
    info.TransformName = 'Sform';
    info.Transform.T(1,1) = params.gen.res(1)*1000;
    info.Transform.T(2,2) = params.gen.res(2)*1000;
    info.Transform.T(3,3) = params.gen.res(3)*1000;
    info.Filemoddate = string(datetime('today','Format','y-MM-d'));

    tmp = strjoin(contrasts,'');

    % Rotate the file if needed
    if rotate
        merged = rot90(merged,-1);
        info.Transform.T(1,1) = params.gen.res(2)*1000;
        info.Transform.T(2,2) = params.gen.res(1)*1000;
        info.ImageSize(1) = params.gen.n(2);
        info.ImageSize(2) = params.gen.n(1);
    end
    
    % Saving NIFTI
    if is2d == 1
        if isempty(b0)
            niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_2d3d_%s%s%s%s.nii',root,folder,scan,tmp,traj,p_dork,r_dork,drift),info)
    %         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_%s_mrreco.nii',root,folder,scan,traj))
        else
            niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_2d3d_%s%s%s%s_b0_%s.nii',root,folder,scan,tmp,traj,p_dork,r_dork,drift,b0),info)
    %         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_b0_%s_%s_mrreco.nii',root,folder,scan,traj,b0))
        end
    else
        if isempty(b0)
            niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_3d_%s%s%s%s.nii',root,folder,scan,tmp,traj,p_dork,r_dork,drift),info)
    %         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_mrreco.nii',root,folder,scan,traj))
        else
            niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_3d_%s%s%s%s_b0_%s.nii',root,folder,scan,tmp,traj,p_dork,r_dork,drift,b0),info)
    %         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_b0_%s_mrreco.nii',root,folder,scan,traj,b0))
        end
    end
    
elseif contains(scan,'cv')
    b0 = '';           % b0='','fessler','romeo','gilad','skope','fessler'
    traj = '';           % 'nom', 'sk'
    is2d = 0;
    p_dork = "";              % Partition DORK "_pDORK"
    r_dork = "";              % Repetition DORK "_rDORK"
    drift = "";               % Drift "_drift"
    
    niftis = dir(fullfile(sprintf('%s/data/%s/raw/nifti/',root,folder),'*nii'));
    niftis = string({niftis.name}');
    idx = contains(niftis,scan);
    niftis = niftis(idx);
   
    vaso = niftiread(sprintf('%s/data/%s/raw/nifti/%s',root,folder,niftis(1)));
    bold = niftiread(sprintf('%s/data/%s/raw/nifti/%s',root,folder,niftis(2)));
    
    info = niftiinfo(sprintf('%s/data/%s/raw/nifti/%s',root,folder,niftis(1)));
    
    
    merged = cat(4,vaso,bold);
    
    merged(:,:,:,1:2:end) = bold;
    merged(:,:,:,2:2:end) = vaso;
    
    info.Filename = sprintf('%s/data/%s/recon/%s_bv_epi.nii',root,folder,scan);
    info.ImageSize(4) = size(merged,4);
    info.raw.dim(5) = size(merged,4);
    
    niftiwrite(merged,sprintf('%s/data/%s/recon/%s_bv_epi.nii',root,folder,scan),info)

end





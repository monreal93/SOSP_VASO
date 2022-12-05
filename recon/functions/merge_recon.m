clear all; clc

cd /mnt/5T3/Alejandro/sosp_vaso/
addpath(genpath('../tools/as'))

folder = '10292022_sv_2';
scan = 'sv_01';
b0 = 'romeo';           % b0='','romeo','gilad','skope'
traj = 'nom';           % 'nom', 'sk'
reps = 1;                % If its only one rep=1, it will use rep2 and won't do time series...
rep_recon = 40;          % Rep to recon if reps =1
is2d = 1;
contrasts = ["b"];  % "v","b","abc"
rotate = false;

% root = '/mnt/ssh/var/www/dabeast/5T3/Alejandro/sosp_vaso';
root = '/mnt/5T3/Alejandro/sosp_vaso';

load(sprintf('%s/data/%s/acq/%s_params.mat',root,folder,scan))


if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
%     contrasts = ["v","b"];
    merged = zeros([params.gen.n reps*2]);
else
%     contrasts = "abc";
    merged = zeros([params.gen.n reps]);
end

k = 1;
for j=1:reps
    if reps == 1
        j = rep_recon;
    end
    if is2d == 1
        for i = 1:params.gen.n(3)
                if  isempty(b0)
                    tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s_mrreco.nii',root,folder,scan,contrasts(1),i,j,traj));
                    merged(:,:,i,k) = tmp;
                    if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                        tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s_mrreco.nii',root,folder,scan,contrasts(2),i,j,traj));
                        merged(:,:,i,k+1) = tmp;
                    end
                else
                    tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s_b0_%s_mrreco.nii',root,folder,scan,contrasts(1),i,j,traj,b0));
                    merged(:,:,i,k) = tmp;
                    if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                        tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_%s_b0_%s_mrreco.nii',root,folder,scan,contrasts(2),i,j,traj,b0));
                        merged(:,:,i,k+1) = tmp;
                    end
                end
        end
    else
        if  isempty(b0)
                tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s_mrreco.nii',root,folder,scan,contrasts(1),j,traj));
                merged(:,:,:,k) = tmp;
                if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                    tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s_mrreco.nii',root,folder,scan,contrasts(2),j,traj));
                    merged(:,:,:,k+1) = tmp;
                end
        else
                tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s_b0_%s_mrreco.nii',root,folder,scan,contrasts(1),j,traj,b0));
                merged(:,:,:,k) = tmp;
                if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
                    tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_%s_b0_%s_mrreco.nii',root,folder,scan,contrasts(2),j,traj,b0));
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

merged = merged.*10e3;

% % Sometimes I need to scale it even more
merged = merged.*1e5;

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
if reps > 1
    info.PixelDimensions(4) = params.gen.volTR;
    if reps>1 && length(contrasts) >1
        info.ImageSize(4) = reps*2;
    else
        info.ImageSize(4) = reps;
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


if is2d == 1
    if isempty(b0)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_2d3d_%s_mrreco.nii',root,folder,scan,tmp,traj),info)
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_%s_mrreco.nii',root,folder,scan,traj))
    else
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_2d3d_b0_%s_%s_mrreco.nii',root,folder,scan,tmp,traj,b0),info)
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_b0_%s_%s_mrreco.nii',root,folder,scan,traj,b0))
    end
else
    if isempty(b0)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_3d_%s_mrreco.nii',root,folder,scan,tmp,traj),info)
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_mrreco.nii',root,folder,scan,traj))
    else
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_%s_3d_%s_b0_%s_mrreco.nii',root,folder,scan,tmp,traj,b0),info)
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_b0_%s_mrreco.nii',root,folder,scan,traj,b0))
    end
end

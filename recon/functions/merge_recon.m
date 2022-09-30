clear all; clc

cd /mnt/5T3/Alejandro/sosp_vaso/
addpath(genpath('../tools/as'))

folder = '09192022_sv';
scan = 'sv_05';
b0 = '';           % b0='','romeo','gilad','skope'
traj = 'nom';           % 'nom', 'sk'
rep = 12;                % If its only one rep=1, it will use rep2 and won't do time series...
is2d = 0;
contrasts = ["b"];  % "v","b","abc"

% root = '/mnt/ssh/var/www/dabeast/5T3/Alejandro/sosp_vaso';
root = '/mnt/5T3/Alejandro/sosp_vaso';

load(sprintf('%s/data/%s/acq/%s_params.mat',root,folder,scan))


if (contains(scan,'sv') || contains(scan,'cv')) && length(contrasts) == 2
%     contrasts = ["v","b"];
    merged = zeros([params.gen.n rep*2]);
else
%     contrasts = "abc";
    merged = zeros([params.gen.n rep]);
end

k = 1;
for j=1:rep
    if rep == 1
        j = j+1;
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

merged = merged.*1e3;

% Adding nifti information
% info = niftiinfo(sprintf('./data/%s/recon/3d/%s_v_rep%i_3d_mrreco.nii',folder,scan,1));
info = niftiinfo(sprintf('./data/%s/recon/3d/%s_v_rep%i_3d_nom_mrreco.nii',folder,'sv_01',2));
info.PixelDimensions = params.gen.res*1000;
info.PixelDimensions(4) = params.gen.volTR;
info.ImageSize = params.gen.n;
info.ImageSize(3) = info.ImageSize(3).*params.gen.kz;
if rep>1
    info.ImageSize(4) = rep*2;
else
    info.ImageSize(4) = rep;
end
info.TimeUnits = 'Second';
info.Datatype = 'single';
info.TransformName = 'Sform';
info.Transform.T(1,1) = params.gen.res(1)*1000;
info.Transform.T(2,2) = params.gen.res(2)*1000;
info.Transform.T(3,3) = params.gen.res(3)*1000;

if is2d == 1
    if isempty(b0)
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_%s_mrreco.nii',root,folder,scan,traj),info)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_%s_mrreco.nii',root,folder,scan,traj))
    else
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_b0_%s_%s_mrreco.nii',root,folder,scan,traj,b0),info)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_b0_%s_%s_mrreco.nii',root,folder,scan,traj,b0))
    end
else
    if isempty(b0)
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_mrreco.nii',root,folder,scan,traj),info)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_mrreco.nii',root,folder,scan,traj))
    else
%         niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_b0_%s_mrreco.nii',root,folder,scan,traj,b0),info)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_%s_b0_%s_mrreco.nii',root,folder,scan,traj,b0))
    end
end

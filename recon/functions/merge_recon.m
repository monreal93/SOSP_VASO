clear all; clc

cd /mnt/5T3/Alejandro/sosp_vaso/
addpath(genpath('../tools/as'))

folder = '09072022_sv_abc';
scan = 'sv_04';
b0 = '';           % b0='','romeo','gilad','skope'
rep = 15;
is2d = 0;

% root = '/mnt/ssh/var/www/dabeast/5T3/Alejandro/sosp_vaso';
root = '/mnt/5T3/Alejandro/sosp_vaso';

load(sprintf('%s/data/%s/acq/%s_params.mat',root,folder,scan))



if contains(scan,'sv') || contains(scan,'cv')
    contrasts = ["v","b"];
    merged = zeros([params.gen.n rep*2]);
else
    contrasts = "abc";
    merged = zeros([params.gen.n rep]);
end

k = 1;
for j=1:rep
    if is2d == 1
        for i = 1:params.gen.n(3)
                if  isempty(b0)
                    tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_mrreco.nii',root,folder,scan,contrasts(1),i,j));
                    merged(:,:,i,k) = tmp;
                    if contains(scan,'sv') || contains(scan,'cv')
                        tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_mrreco.nii',root,folder,scan,contrasts(2),i,j));
                        merged(:,:,i,k+1) = tmp;
                    end
                else
                    tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_b0_%s_mrreco.nii',root,folder,scan,contrasts(1),i,j,b0));
                    merged(:,:,i,k) = tmp;
                    if contains(scan,'sv') || contains(scan,'cv')
                        tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_%s_sl%i_rep%i_2d_b0_%s_mrreco.nii',root,folder,scan,contrasts(2),i,j,b0));
                        merged(:,:,i,k+1) = tmp;
                    end
                end
        end
    else
        if  isempty(b0)
                tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_mrreco.nii',root,folder,scan,contrasts(1),j));
                merged(:,:,:,k) = tmp;
                if contains(scan,'sv') || contains(scan,'cv')
                    tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_mrreco.nii',root,folder,scan,contrasts(2),j));
                    merged(:,:,:,k+1) = tmp;
                end
        else
                tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_b0_%s_mrreco.nii',root,folder,scan,contrasts(1),j,b0));
                merged(:,:,:,k) = tmp;
                if contains(scan,'sv') || contains(scan,'cv')
                    tmp = niftiread(sprintf('%s/data/%s/recon/3d/%s_%s_rep%i_3d_b0_%s_mrreco.nii',root,folder,scan,contrasts(2),j,b0));
                    merged(:,:,:,k+1) = tmp;
                end
        end
    end
if contains(scan,'sv') || contains(scan,'cv')
    k = k+2;
else
    k = k+1;
end
end

if is2d == 1
    if isempty(b0)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_mrreco.nii',root,folder,scan))
    else
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_2d3d_b0_%s_mrreco.nii',root,folder,scan,b0))
    end
else
    if isempty(b0)
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_mrreco.nii',root,folder,scan))
    else
        niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_3d_b0_%s_mrreco.nii',root,folder,scan,b0))
    end
end


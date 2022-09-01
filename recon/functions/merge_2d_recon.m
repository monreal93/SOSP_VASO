clear all; clc

folder = '08312022_sv';
scan = 'sv_05';
b0 = 'romeo';           % b0='','romeo','gilad','skope'
% contrasts = ['v','b'];      % 'v','b'
rep = 120;

root = '/mnt/ssh/var/www/dabeast/5T3/Alejandro/sosp_vaso';

load(sprintf('%s/data/%s/acq/%s_params.mat',root,folder,scan))

merged = zeros([params.gen.n rep*2]);

k = 1;
for j=1:rep
    for i = 1:params.gen.n(3)
            if  isempty(b0)
                tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_v_sl%i_rep2_2d_mrreco.nii',root,folder,scan,i));
            else
                tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_v_sl%i_rep%i_2d_b0_%s_mrreco.nii',root,folder,scan,i,j,b0));
                merged(:,:,i,k) = tmp;
                tmp = niftiread(sprintf('%s/data/%s/recon/2d/%s_b_sl%i_rep%i_2d_b0_%s_mrreco.nii',root,folder,scan,i,j,b0));
                merged(:,:,i,k+1) = tmp;
            end
    end
k = k+2;
end

niftiwrite(single(merged),sprintf('%s/data/%s/recon/%s_b0_%s_mrreco.nii',root,folder,scan,b0))



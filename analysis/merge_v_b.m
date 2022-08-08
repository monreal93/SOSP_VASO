
% cd /mnt/ssh/var/www/dabeast/5T3/Alejandro/sosp_vaso/

cd /mnt/5T3/Alejandro/sosp_vaso/

scan = 'sv_01';
folder = '08052022_sv_abc';
rep = 120;

% tmp = sprintf('./data/%s/recon/*%s*',folder,scan);
% files = dir(tmp);

for i=1:rep

%     tmp = getfield(files,{i},'name');
    tmp = sprintf('./data/%s/recon/%s_v_rep%i_3d_mrreco.nii',folder,scan,i);
    vaso(:,:,:,i) = niftiread(tmp);

    tmp = sprintf('./data/%s/recon/%s_b_rep%i_3d_mrreco.nii',folder,scan,i);
    bold(:,:,:,i) = niftiread(tmp);
   
end

v_b = zeros(size(vaso,1),size(vaso,2),size(vaso,3)*2,size(vaso,4));
v_b(:,:,1:2:end,:) = vaso;
v_b(:,:,2:2:end,:) = bold;

niftiwrite(v_b,sprintf('./data/%s/analysis/%s_v_b.nii',folder,scan))


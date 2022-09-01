
% cd /mnt/ssh/var/www/dabeast/5T3/Alejandro/sosp_vaso/

clear

cd /mnt/5T3/Alejandro/sosp_vaso/

scan = 'sv_01';
folder = '08312022_sv';
rep = 120;  % 120

% tmp = sprintf('./data/%s/recon/*%s*',folder,scan);
% files = dir(tmp);

load(sprintf('./data/%s/acq/%s_params',folder,scan))

% v_b = zeros(size(vaso,1),size(vaso,2),size(vaso,3),size(vaso,4)*2);
j = 1;
for i=1:rep

%     % This way I am doing V-B
%     tmp = sprintf('./data/%s/recon/%s_v_rep%i_3d_mrreco.nii',folder,scan,i);
%     v = niftiread(tmp);
%     v_b(:,:,:,j) = v;
%     
%     tmp = sprintf('./data/%s/recon/%s_b_rep%i_3d_mrreco.nii',folder,scan,i);
%     b = niftiread(tmp);
%     v_b(:,:,:,j+1) = b;
%     
%     j = j+2;

    % This way I am doing B-V
    tmp = sprintf('./data/%s/recon/3d/%s_v_rep%i_3d_mrreco.nii',folder,scan,i);
    v = niftiread(tmp);
    v_b(:,:,:,j+1) = v;
    
    tmp = sprintf('./data/%s/recon/3d/%s_b_rep%i_3d_mrreco.nii',folder,scan,i);
    b = niftiread(tmp);
    v_b(:,:,:,j) = b;
    
    j = j+2;
end

% v_b = zeros(size(vaso,1),size(vaso,2),size(vaso,3)*2,size(vaso,4));
% v_b(:,:,1:2:end,:) = vaso;
% v_b(:,:,2:2:end,:) = bold;

% Scaling since values are too small
v_b = v_b.*1e4;

% Adding nifti information
info = niftiinfo(sprintf('./data/%s/recon/%s_v_rep%i_3d_mrreco.nii',folder,scan,1));
info.PixelDimensions = params.gen.res*1000;
info.PixelDimensions(4) = params.gen.volTR;
info.ImageSize = params.gen.n;
info.ImageSize(3) = info.ImageSize(3).*params.gen.kz;
info.ImageSize(4) = rep*2;
info.TimeUnits = 'Second';
info.Datatype = 'single';
info.TransformName = 'Sform';
info.Transform.T(1,1) = params.gen.res(1)*1000;
info.Transform.T(2,2) = params.gen.res(2)*1000;
info.Transform.T(3,3) = params.gen.res(3)*1000;

niftiwrite(single(v_b),sprintf('./data/%s/analysis/%s_v_b.nii',folder,scan),info)


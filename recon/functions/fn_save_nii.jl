function fn_save_nii(volume,save_name)

    vol = NIVolume(abs.(volume));
    path = string("/usr/share/sosp_vaso/data/tmp/",save_name,".nii");
    niwrite(path,vol);
end
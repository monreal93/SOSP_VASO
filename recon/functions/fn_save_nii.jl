function fn_save_nii(volume,save_name)

    vol = NIVolume(real.(volume));
    path = string("/usr/share/sosp_vaso/data/tmp/",save_name,".nii");
    niwrite(path,vol);
end
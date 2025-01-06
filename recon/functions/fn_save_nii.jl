function fn_save_nii(volume,save_name)

    vol = NIVolume(real.(volume));
    path = string("/neurodesktop-storage/5T4/Alejandro/sosp_vaso/data/tmp/",save_name,".nii");
    niwrite(path,vol);
end
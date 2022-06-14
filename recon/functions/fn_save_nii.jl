macro Name(arg)
    string(arg)
end

function fn_save_nii(volume)

    vol = NIVolume(abs.(volume));
    path = string("/usr/share/Alejandro/tmp/data/",@Name(volume),".nii");
    niwrite(path,vol)
end
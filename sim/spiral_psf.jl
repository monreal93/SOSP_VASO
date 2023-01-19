using MAT, MRIReco, NIfTI
using infiltrator, MIRTjim

include("../recon/functions/fn_save_nii.jl")

scans = ["sv_12"]#,"sv_02","sv_03","sv_04","sv_05","sv_06","sv_07"]; #,"sv_02","sv_03","sv_04","sv_05"];
directory = "11252022_psf"        # directory where the data is stored
obj = 1;                          # Simulate 0=point, 1=brain

mtx_s = (238,240,24)


if obj == 1
    img = niread(string("/usr/share/sosp_vaso/data/",directory,"/tmp/gre_238_240_24_ech1.nii"))
    img = img.raw
else
    img = zeros(mtx_s)
    img[Int(mtx_s[1]/2),Int(mtx_s[2]/2),Int(mtx_s[3]/2)] = 1;
end

@infiltrate
for i = 1:Int(length(scans))

    params =  matread(string("/usr/share/sosp_vaso/data/",directory,"/acq/",scans[i],"_params.mat"))
    params = params["params"]
    traj = matread(string("/usr/share/sosp_vaso/data/",directory,"/acq/",scans[i],"_ks_traj_nom.mat"))
    traj = traj["ks_traj"]

    # mtx_s = (Tuple{Int,Int,Int}(params["gen"]["n"]))

    traj = [traj["kx"][:] traj["ky"][:] traj["kz"][:]]
    traj = permutedims(traj,[2,1])

    # Normalizing the traj to -0.5 to 0.5
    kx_fov = maximum([abs.(minimum(traj[1,:])) abs.(maximum(traj[1,:]))])*2
    ky_fov = maximum([abs.(minimum(traj[2,:])) abs.(maximum(traj[2,:]))])*2
    kz_fov = maximum([abs.(minimum(traj[3,:])) abs.(maximum(traj[3,:]))])*2

    traj[1,:] = traj[1,:]./kx_fov 
    traj[2,:] = traj[2,:]./ky_fov 
    traj[3,:] = traj[3,:]./kz_fov 

    op = NFFTOp(mtx_s,traj)

    img_ft = op*vec(img)

    img_new = adjoint(op)*img_ft
    img_new = reshape(img_new,mtx_s)
    img_new = abs.(img_new)
    img_new = convert(Array{Float32}, img_new)

    if obj == 1
        matwrite(string("/usr/share/sosp_vaso/data/",directory,"/tmp/",scans[i],"_sim_recon.mat"),Dict("recon" => img_new))
        vol = NIVolume(abs.(img_new))
        path = string("/usr/share/sosp_vaso/data/",directory,"/tmp/",scans[i],"_sim_recon.nii")
        niwrite(path,vol)
    else
        matwrite(string("/usr/share/sosp_vaso/data/",directory,"/acq/",scans[i],"_psf.mat"),Dict("psf" => img_new))
        vol = NIVolume(abs.(img_new))
        path = string("/usr/share/sosp_vaso/data/",directory,"/tmp/",scans[i],"_psf.nii")
        niwrite(path,vol)
    end

    # fn_save_nii(img_new,string("psf_",scans[i]))
end


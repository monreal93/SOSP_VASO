# using Infiltrator
using SphericalHarmonicExpansions, Distributions
# using MIRTjim: jim
# using MAT
# using ImageTransformations

function getDynamicOffResonanceCalibrationMatrix(params::Dict{Symbol,Any},b0_init,recon_b0; calib_steps=100)
    
    n_channels = size(recon_b0,4)
    n_times = calib_steps

    recon_b0 = recon_b0[:,:,:,:,1]
    recon_b0 = imresize(recon_b0,size(b0_init)[1:3])

    # b0_init = imresize(b0_init,size(recon_b0)[1:3])

    mtx_s = size(recon_b0)[1:3]
    fov = params_pulseq["gen"]["fov"]
    l = 2
    lmax = (l+1).^2

    # b0_init = matread("../data/04062023_sv_abc/acq/b0_238_236_24.mat");
    # b0_init = b0_init["b0"]
    # b0_init = b0

    # Convert B0_init from rad/s to Hz
    b0_init = imag.(b0_init)./(2*pi)

    sh_basis = zeros(mtx_s[1],mtx_s[2],mtx_s[3],lmax)

    # if params[:recon_b0ri] == 1
    #     recon_b0 = matread(string("../",params[:directory],"acq/recon_b0_", params[:scan],".mat"))
    # else
    #     recon_b0 = matread(string("../",params[:directory],"acq/recon_b0_",params[:scan][1:2],"_01.mat"))
    # end
    # recon_b0 = recon_b0["fieldmap"][:,:,:,2,:]

    # recon_b0 = imresize(recon_b0,(mtx_s[1],mtx_s[2],mtx_s[3],params[:nCoils]))
    # recon_b0 = recon_b0[:,:,:,:,1]
    # recon_b0 = recon_b0.*1e-4

    # Getting shperical harmonics basis functions
    @polyvar x y z

    l = 0:l
    m = -l[end]:l[end]

    # Creating cartesian coordinates
    (xx,yy,zz) = (LinRange.(-1,1,[mtx_s[1],mtx_s[2],mtx_s[3]]).*vec(fov)./2)

    tmp = 1
    for  il=0:l[end]
        for im=-il:il
                f = rlylm(il,im,x,y,z)
                g = @fastfunc f
                for ix=eachindex(xx), iy=eachindex(yy), iz=eachindex(zz)
                    sh_basis[ix, iy, iz,tmp] = Base.invokelatest(g,xx[ix], yy[iy], zz[iz])
                end

                tmp = tmp+1
                @info string("l=",il,", m=",im)
        end
    end

    # Now trying to "decompose" the initial B0 map into its shperical harmonics
    # b0 = reshape(b0,:,lmax)
    # b0_init = vec(imag.(b0_init))

    sh_basis = reshape(sh_basis,:,lmax)

    b = sh_basis \ vec(b0_init)

    # del_b0 is the lth-order aproximation of the initial b0 map
    ΔB0 = sh_basis * b
    ΔB0 = reshape(ΔB0,mtx_s[1],mtx_s[2],mtx_s[3])

    b_calib = zeros(lmax,calib_steps)
    # Random delta B0 coefficients, trying to get a realistic order of magnitude using ΔB0 as reference +/-10%
    for i=1:calib_steps    
        for j=1:lmax    
            # d = Uniform(-abs(b[j]).*0.1,abs(b[j]).*0.1)
            d = Uniform(0,abs(b[j]).*0.1)
            # d = Uniform(0,1)
            b_calib[j,i] = rand(d)
        end
    end

    fid = Matrix{ComplexF32}(undef,n_channels,n_times)
    calib = Matrix{ComplexF32}(undef,n_channels,lmax)

    @floop for i=1:calib_steps

        sh_tmp = sh_basis * b_calib[:,i]

        sh_tmp = reshape(sh_tmp,mtx_s[1],mtx_s[2],mtx_s[3])

        img_sim = recon_b0.*exp.(-im*2*pi*2e-3.*sh_tmp).*exp.(im*2*pi*2e-3.*b0_init)

        fid[:,i] = dropdims(sum(img_sim,dims=(1,2,3)),dims=(1,2,3))

    end

    calib = fid / b_calib

    @info("Stop... Dyn B0 corr...")
    @infiltrate

    # # Now generating the calibration matrix A as im Wallace paper
    # s_vec = reshape(recon_b0,:,params[:nCoils])
    # b0_vec = reshape(sh_basis,:,lmax)

    # A = Array{ComplexF32, 2}(undef, params[:nCoils],lmax)
    # B = Array{ComplexF32, 2}(undef, mtx_s[1]*mtx_s[2]*mtx_s[3],lmax)

    # # Generate A (calibration) matrix
    # for i=1:lmax, j=1:params[:nCoils]
    #     A[j,i] = s_vec[:,j]'*exp.(im.*(42.58e6).*params[:b0_te][1]*1e-3.*b0_vec[:,i]*b[i])
    #     # A[j,i] = s_vec[:,j]'*exp.(im.*params[:b0_te][1].*b0_vec[:,i]*b[i])
    # end

    # # Getting sh decomposition of ΔB0
    # for i=1:lmax
    #     tmp = b0_vec[:,i]*b[i]
    #     # tmp = reshape(tmp,238,236,24)
    #     B[:,i] =  tmp
    # end

    return calib, sh_basis
end

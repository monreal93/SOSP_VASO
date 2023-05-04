# using Infiltrator
using SphericalHarmonicExpansions
# using MIRTjim: jim
# using MAT
# using ImageTransformations

function getCalibrationMatrix(params::Dict{Symbol,Any},b0_init)
    mtx_s = [Int64(params[:mtx_s][1]),Int64(params[:mtx_s][2]),Int64(params[:mtx_s][3])]
    fov = params[:fov]
    l = 3
    lmax = (l+1).^2


    # b0_init = matread("../data/04062023_sv_abc/acq/b0_238_236_24.mat");
    # b0_init = b0_init["b0"]
    # b0_init = b0

    # Convert B0_init from rad/s to Hz
    b0_init = imag.(b0_init)./(2*pi)

    sh_basis = zeros(mtx_s[1],mtx_s[2],mtx_s[3],lmax)

    if params[:fmri] == 1
        fm = matread(string("../",params[:directory],"acq/fm_", params[:scan],".mat"))
    else
        fm = matread(string("../",params[:directory],"acq/fm_",params[:scan][1:2],"_01.mat"))
    end
    
    fm = fm["fieldmap"][:,:,:,3,:]
    fm = imresize(fm,(mtx_s[1],mtx_s[2],mtx_s[3],params[:nCoils]))
    # fm = abs.(fm)


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

    # Now generating the calibration matrix A as im Wallace paper
    s_vec = reshape(fm,:,params[:nCoils])
    b0_vec = reshape(sh_basis,:,lmax)

    A = Array{ComplexF32, 2}(undef, params[:nCoils],lmax)
    B = Array{ComplexF32, 2}(undef, mtx_s[1]*mtx_s[2]*mtx_s[3],lmax)

    # Generate A (calibration) matrix
    for i=1:lmax, j=1:params[:nCoils]
        A[j,i] = s_vec[:,j]'*exp.(im.*(42.58e6).*params[:b0_te][1]*1e-3.*b0_vec[:,i]*b[i])
        # A[j,i] = s_vec[:,j]'*exp.(im.*params[:b0_te][1].*b0_vec[:,i]*b[i])
    end

    # Getting sh decomposition of ΔB0
    for i=1:lmax
        tmp = b0_vec[:,i]*b[i]
        # tmp = reshape(tmp,238,236,24)
        B[:,i] =  tmp
    end
    
    jim(ΔB0, "ΔB0"; color=:jet)
    @infiltrate

    return A, B, sh_basis, ΔB0, b
end

using SphericalHarmonicExpansions

"""
CalculateSphericalHarmonics(params::Dict{Symbol,Any},SensitivityMap,OffResonanceMap,ks_traj; l::Int=3)

Creates Spherical harmonic basis functions
"""
function CalculateSphericalHarmonics(params::Dict{Symbol,Any},SensitivityMap,OffResonanceMap,ks_traj; l::Int=3)
    reconSize = [Int64(params[:reconSize][1]),Int64(params[:reconSize][2]),Int64(params[:reconSize][3])]
    fov = params[:FieldOfView]
    lmax = (l+1).^2

    SphericalHarmonicsBasis = zeros(reconSize[1],reconSize[2],reconSize[3],lmax)


    fm = imresize(fm,(reconSize[1],reconSize[2],reconSize[3],params[:nCoils]))
    fm = fm[:,:,:,1,:]
    # fm = fm.*1e-4

    # Getting shperical harmonics basis functions
    @polyvar x y z

    l = 0:l
    m = -l[end]:l[end]

    # Creating cartesian coordinates
    (xx,yy,zz) = (LinRange.(-1,1,[reconSize[1],reconSize[2],reconSize[3]]).*vec(fov)./2)

    tmp = 1
    for  il=0:l[end]
        for im=-il:il
                f = rlylm(il,im,x,y,z)
                g = @fastfunc f
                for ix=eachindex(xx), iy=eachindex(yy), iz=eachindex(zz)
                    SphericalHarmonicsBasis[ix, iy, iz,tmp] = Base.invokelatest(g,xx[ix], yy[iy], zz[iz])
                end

                tmp = tmp+1
                @info string("l=",il,", m=",im)
        end
    end

    SphericalHarmonicsBasis = reshape(SphericalHarmonicsBasis,:,lmax)

    return SphericalHarmonicsBasis
end


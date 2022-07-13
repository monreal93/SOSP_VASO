
function fn_ismrmd(params_sv::Dict{Symbol,Any})

    

    o_file = string(params_sv[:path],"ismrmd/",params_sv[:scan],".mrd");

    # First we delete the ismrmd file if it already exists
    cmd = string("rm -f ",o_file);
    cmd = `sh -c $cmd`;
    run(cmd)

    # siemens to ismrmd convertion
    cmd = string("/siemens_to_ismrmrd/build/siemens_to_ismrmrd -f ",params_sv[:path],"raw/twix/*",
        params_sv[:scan],"* -o ",o_file);
    cmd = `sh -c $cmd`;
    run(cmd)

   
    # Reading and reshaping the just created ismrmd file
    file = ISMRMRDFile(o_file);
    rawData = RawAcquisitionData(file);

    # First I concatenate all the kdata
    tmp = rawData.profiles[1].data;
    for i in 2:first(collect(size(rawData.profiles))) 
        tmp = vcat(tmp,rawData.profiles[i].data)
    end



    @infiltrate

end
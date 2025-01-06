function gz_blips = prepare_gz_blips(params)
lims = params.gen.lims;    

% if params.gen.kz_caipi > 1    
%     i_tmp = params.gen.n(3)/params.gen.seg;
% else
    i_tmp = params.gen.n(3);
% end
% Kz Blips
tmp = 1;
for i=1:i_tmp
    for j=1:params.gen.seg
        if mod(params.gen.n(3),2) == 1
            area = -(params.gen.del_k(3)*(ceil(params.gen.n(3)/2)+1))+(params.gen.del_k(3)/params.gen.kz_caipi*(tmp-1))+(params.gen.del_k(3)/2);
        else
            area = -(params.gen.del_k(3)*(ceil(params.gen.n(3)/2)))+(params.gen.del_k(3)/params.gen.kz_caipi*(tmp-1));
        end
        if area==0
            area = 1e-10;
        end
        dur = ceil(2*sqrt(area/lims.maxSlew)/10e-6)*10e-6;
        if i == 1
            gz_blips(i,j) = mr.makeTrapezoid('z',lims,'Area',area);
        else
            gz_blips(i,j) = mr.makeTrapezoid('z',lims,'Area',area,'Duration',mr.calcDuration(gz_blips(1)));
        end
        if params.gen.kz_caipi > 1
            tmp = tmp+1;
        end
    end
    if params.gen.kz_caipi == 1
        tmp = tmp+1;
    end
end

% ToDo: Fix this....
% Reshuffling blips if center-out
if params.gen.kz_enc == 1
    tmp = [];
    j = floor(floor(params.gen.n(3)/params.gen.kz)/2);
    for i=1:length(gz_blips)
        tmp = [tmp gz_blips(j)];
        if mod(i,2) == 0
            j = j-i;
        else
            j = j+i;
        end
    end
    gz_blips = tmp;
end

end
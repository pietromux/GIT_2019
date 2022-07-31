fid=fopen(totalmagfilename,'r');
tline = fgets(fid);
awind = 1;
qfind = 1;
adind = 1;
lat = zeros(3,1);
while ischar(tline)
    tline = fgets(fid);
    if strfind(tline,'AW')
        tline(1:2) = [];
        aw = str2num(tline);
        lat(1,aw(3)+awind : awind+aw(3)+aw(2)-1) = aw(1);
        awind = awind+aw(3)+aw(2);
    end
    if strfind(tline,'QF')
        tline(1:2) = [];
        qf = str2num(tline);
        lat(2,qf(3)+qfind:qfind+qf(3)+qf(2)-1) = qf(1);
        qfind = qfind+qf(3)+qf(2);
    end
    if strfind(tline,'AD')
        tline(1:2) = [];
        ad = str2num(tline);
        lat(3,adind+ad(3):adind+ad(3)+ad(2)-1) = ad(1);
        adind= adind+ad(3)+ad(2);
    end
end
fclose(fid);
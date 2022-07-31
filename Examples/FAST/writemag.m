quadgradient = 160;
lamu = 0.032;
nwig = 1;

totalmagfilename = 'mg_und.mag.in';
readin_totalmagfile;
lat(2,lat(2,:)==283.4000)=quadgradient;
lat(2,lat(2,:)==-283.4000)=-quadgradient;

fname=['mg_undnew.mag.in'];
fid=fopen(fname,'w');
fprintf(fid,'%s\n','? VERSION=1.0');
fprintf(fid,'%s %f\n','? UNITLENGTH=',lamu);
for j = 1:size(lat,2)
    fprintf(fid,'%s %f %.f %.f\n', 'AW',lat(1,j),nwig,0);
    fprintf(fid,'%s %.3f %.f %.f\n','QF',lat(2,j),nwig,0);
    fprintf(fid,'%s %.3f %.f %.f\n','AD',lat(3,j),nwig,0);
end

fclose(fid);

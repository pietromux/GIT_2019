function fid = writemag(z,x,sep,aw0,B0,lamu,nwig)
      
fname=['C:\cygwin\home\pietro\TESS3.0\','mg','.mag.in'];

fid=fopen(fname,'w');

undd=0;
quadd=0;

fprintf(fid,'%s\n','? VERSION=1.0');
fprintf(fid,'%s %f\n','? UNITLENGTH=',lamu);

xnum=z2x(z,x,sep,lamu);

switch xnum
    
    case 5
    fprintf(fid,'%s %f %.f %.f\n', 'AW',aw0, nwig, undd);
    fprintf(fid,'%s %.3f %.f %.f\n','QF',B0,nwig,quadd);
    case 2
    fprintf(fid,'%s %f %.f %.f\n', 'AW', aw0, nwig, undd);
    case 6
    fprintf(fid,'%s %f %.f %.f\n', 'AW', aw0, nwig, undd);
    fprintf(fid,'%s %.3f %.f %.f\n','QF',-B0,nwig,quadd);
    case 1
    fprintf(fid,'%s %f %.f %.f\n', 'AW', 0 ,nwig, undd);

end
 fprintf(fid,'%s %f %.f %.f\n', 'AW', 0 ,0, 0);
fclose(fid);
end
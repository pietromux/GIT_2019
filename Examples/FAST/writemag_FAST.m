
      
fname=['C:\cygwin\home\pietro\TESS3.0\','mg','.mag.in'];

fid=fopen(fname,'w');

undd=0;
quadd=0;

fprintf(fid,'%s\n','? VERSION=1.0');
fprintf(fid,'%s %f\n','? UNITLENGTH=',lamu);

    fprintf(fid,'%s %f %.f %.f\n', 'AW',aw0, nwig, undd);
    fprintf(fid,'%s %.3f %.f %.f\n','QF',B0,nwig,quadd);

    fclose(fid);

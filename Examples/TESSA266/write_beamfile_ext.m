beam= importdata('PTvary007.twiss4.beam',' ',2);
filename = strcat('constantenergy.beam');
fid = fopen(filename,'wt');
fprintf(fid, ' ? VERSION 1.0 \n');
onlylps = 0;
if(onlylps)
fprintf(fid, ' ? COLUMNS ZPOS CURPEAK GAMMA0 DELGAM \n');
current = smooth(beam.data(:,6),20);
for j = 1:git_par.nslice 
     zpos_beam(j)=git_par.zsep*git_par.lambda*(j-1);
     fprintf(fid, ' %E %E  %E %E\n', zpos_beam(j), current(j), beam.data(j,2), beam.data(j,3));
end
else
fprintf(fid, ' ? COLUMNS ZPOS CURPEAK GAMMA0 DELGAM EMITX EMITY \n');
current = smooth(beam.data(:,6),20);
for j = 1:git_par.nslice 
     zpos_beam(j)=git_par.zsep*git_par.lambda*(j-1);
     fprintf(fid, ' %E %E %E %E %E %E\n', zpos_beam(j), current(j), beam.data(j,2), beam.data(j,3), beam.data(j,4), beam.data(j,5));
end
end 
 

% fprintf(fid, ' ? COLUMNS ZPOS CURPEAK \n');
% for j = 1:git_par.nslice 
%          zpos_beam(j)=git_par.zsep*git_par.lambda*(j-1);
%         fprintf(fid, ' %E %E\n', zpos_beam(j), current(j));
% end
st = fclose('all');

slicenum_rad=1:git_par.nslice;
xlamds=git_par.lambda;
zpos_rad=slicenum_rad*xlamds;
p_rad = git_par.prad0*exp(-((slicenum_rad-radshift)*xlamds).^2/2/(git_par.pulselength^2));

filename = strcat(git_par.dirname,'mg.rad.in');
fid = fopen(filename,'wt');
fprintf(fid, ' ? VERSION 1.0 \n');
fprintf(fid, ' ? COLUMNS ZPOS PRAD0\n');

for jz = 1:length(slicenum_rad)
    fprintf(fid,' %.5e %.5e\n',zpos_rad(jz), p_rad(jz));
end
fclose(fid);

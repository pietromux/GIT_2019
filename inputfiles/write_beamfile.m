%writes beam file in Gaussian distribution, requires curpeakpos (position of current in units
%of #slices) and curlen (FWHM in units of meters.)
clearvars zpos_beam current_norm current_beam
filename = strcat(git_par.dirname,'mg.beam.in');
fid = fopen(filename,'wt');
fprintf(fid, ' ? VERSION 1.0 \n');
fprintf(fid, ' ? COLUMNS ZPOS CURPEAK\n');

if beamtype==0 %%rectangular distribution and get new current length with given charge and current peak
    curlen_slices =git_par.curlen/git_par.lambda;
    current = zeros(1,git_par.nslice);
    current(1,round(git_par.nslice/2-curlen_slices/2):round(git_par.nslice/2+curlen_slices/2))=git_par.curpeak;

elseif beamtype==1 %%Gaussian distribution, (new charge) use curlen and curpeak as is to get gaussian current
    for j = 1:git_par.nslice
   current(j) = git_par.curpeak*exp(-(((j-1)-curpeakpos)*git_par.lambda).^2/2/(git_par.curlen)^2);
    end
    
    git_par.charge=git_par.curpeak/sqrt(1/(2*pi*(git_par.curlen/3e8)^2));
elseif beamtype==2 %%get new current peak with given charge and current length
    for j = 1:git_par.nslice
        current_norm(j)=exp(-(((j-1)-curpeakpos)*git_par.lambda).^2/2/(git_par.curlen)^2);
    end
       newcurpeak=git_par.charge*sqrt(1/(2*pi*(git_par.curlen/3e8)^2));
       git_par.curpeak=newcurpeak;
       current=newcurpeak*current_norm; 

elseif beamtype==3 %%(new charge) use curlen and curpeak as is to get rectangular current
    curlen_slices = git_par.curlen/git_par.lambda;
    current = zeros(1,git_par.nslice);
    current(1,round(git_par.nslice/2-curlen_slices/2):round(git_par.nslice/2+curlen_slices/2))=git_par.curpeak;
    git_par.charge=git_par.curlen/3e8*git_par.curpeak;
    
elseif beamtype==4 %%get new current peak with given charge and current length
    curlen_slices = git_par.curlen/git_par.lambda;
    git_par.curpeak = git_par.charge/(git_par.curlen/3e8);
    current = zeros(1,git_par.nslice);
    current(1,round(git_par.nslice/2-curlen_slices/2):round(git_par.nslice/2+curlen_slices/2))=git_par.curpeak;
end

for j = 1:git_par.nslice 
    zpos_beam(j)=git_par.zsep*git_par.lambda*(j-1);
    fprintf(fid, ' %E %E\n', zpos_beam(j), current(j));
end
st = fclose('all');

% 
% figure(101)
% plot(zpos_beam/git_par.lambda,current)
% xlabel('slicenum')
% ylabel('current(A)')
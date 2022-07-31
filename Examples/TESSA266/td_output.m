zpos=[1:nslice]*git_par.lambda;timearray=[1:nslice]*xlamds/3e8;
Eradf = trapz(timearray,powereff(:,end))
Erad0 = trapz(timearray,powereff(:,1));
Ebeam0 = mean(energy(:,1))*0.511e6*trapz(timearray,current);
Ebeamf = mean(energy(:,end))*0.511e6*trapz(timearray,current);
Eradf/Ebeam0
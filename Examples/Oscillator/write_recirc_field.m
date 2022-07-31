%corrects beam transmission for the time-independent scheme
%For time dependent, edit slippage and filter the radiation according to
%cavity detuning and Q-factor of filter

%% Load field

outputfieldfile = strcat(git_par.dirname,'mg.out.dfl');
Cfieldfile = fopen(outputfieldfile);
Cfield = fread(Cfieldfile,2*git_par.ncar*git_par.ncar*niter,'double');
fclose(Cfieldfile);

%filter with dispersion
Q=git_par.Q;
omega_m=git_par.nslice/2;
for jfreq = 1:git_par.nslice
filter3(jfreq) = 1i*jfreq/Q / (omega_m^2-jfreq^2+1i*jfreq/Q); 
end

%% Transmission
for islice = 1:niter
    iindex = git_par.ncar*git_par.ncar*(islice-1)*2;
    recfield = Cfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1);
    imcfield = Cfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar);
    compfield=recfield+1i*imcfield;
    
    newrecfield=recfield*sqrt(git_par.transmission);
    newimcfield=imcfield*sqrt(git_par.transmission);
    newcompfield=newrecfield+1i*newimcfield;

    newCfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1)=newrecfield;
    newCfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar)=newimcfield;
end

Cfield = newCfield;

if(git_par.itdp)
   
%% Cavity detuning    
newCfield=Cfield*0;
if (git_par.cavitydetuning>0)
    for islice=git_par.cavitydetuning:niter
    iindex = git_par.ncar*git_par.ncar*(islice-1)*2;    
    newiindex=git_par.ncar*git_par.ncar*(islice-git_par.cavitydetuning)*2;  
    newCfield(newiindex+1:2:newiindex+2*git_par.ncar*git_par.ncar-1)= Cfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1);
    newCfield(newiindex+2:2:newiindex+2*git_par.ncar*git_par.ncar)=  Cfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar);
    end
else
    for islice=1:niter+git_par.cavitydetuning
    iindex = git_par.ncar*git_par.ncar*(islice-1)*2;    
    newiindex=git_par.ncar*git_par.ncar*(islice-git_par.cavitydetuning)*2;  
    newCfield(newiindex+1:2:newiindex+2*git_par.ncar*git_par.ncar-1) = Cfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1);
    newCfield(newiindex+2:2:newiindex+2*git_par.ncar*git_par.ncar) =  Cfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar);
    end
end

Cfield = newCfield;

%% Spectral Filtering
 clear recfield imcfield compfield filterfield
for ix = 1:git_par.ncar
    for iy = 1:git_par.ncar
        for islice=1:niter
        iindex = git_par.ncar*git_par.ncar*(islice-1)*2;    
        recfield(islice) = Cfield(iindex+(iy-1)*git_par.ncar*2+ix*2-1);
        imcfield(islice) = Cfield(iindex+(iy-1)*git_par.ncar*2+ix*2);
        compfield(islice) = recfield(islice)-1i*imcfield(islice);
        end
        filterfield = ifft(ifftshift(fftshift(fft(compfield) ).*filter3)); 
%        filterfield = ifft(ifftshift(fftshift(fft(compfield) ) ) );

        newrecfield = real(filterfield);
        newimcfield = imag(filterfield);
        for islice=1:niter
        iindex = git_par.ncar*git_par.ncar*(islice-1)*2;    
        newCfield(iindex+(iy-1)*git_par.ncar*2+ix*2-1) = newrecfield(islice); 
        newCfield(iindex+(iy-1)*git_par.ncar*2+ix*2) = -newimcfield(islice);
        end        
    end
end
end        
%%
outputfieldfile = strcat(git_par.dirname,'mg.out.dfl');
Cfieldfile=fopen(outputfieldfile,'w+');
Ccount = fwrite(Cfieldfile,newCfield,'double');
fclose(Cfieldfile);

% figure(19)
% clf;
% plot(abs(compfield));
% hold on
% plot(abs(newcompfield));
% figure(20)
% clf;
% subplot(2,1,1)
% contourf(abs(reshape(compfield,ncar,ncar)))
% hold on
% subplot(2,1,2)
% contourf(abs(reshape(newcompfield,ncar,ncar)))


%% Open particle and field file for tapering optimization
outputpartfile = strcat(dirname,'mg.out.dpa');

Cpoint = fopen(outputpartfile,'r+');
Ctot = fread(Cpoint,6*npart*niter,'double');
fclose(Cpoint);
bc = 0.12e6;
offset_con= git_par.delgamma/git_par.gamma0/2/bc/git_par.rxbeam
for islice = 1:niter
    for ip =1:npart
        gamma(ip) = Ctot(ip+(islice-1)*6*npart); 
        gamma_avg = mean(gamma);
    end 
    for ip = 1:npart
        rel_gamma(ip) = abs(gamma(ip)-min(gamma))/gamma_avg;
        phiabs(ip) = Ctot(ip+npart+(islice-1)*6*npart);
        phiabs_fin(ip) = phiabs(ip) + 2*pi*R56buncher/lambda*(gamma(ip)-gamma_avg)/gamma_avg+phaseshift;        
        xpos(ip) = Ctot(ip+2*npart+(islice-1)*6*npart)/10+offset_con;
        xpos(ip) = xpos(ip)+sign(xpos(ip))*sqrt(rel_gamma(ip)/bc);
        xprime(ip) = Ctot(ip+4*npart+(islice-1)*6*npart);
        ypos(ip) = Ctot(ip+3*npart+(islice-1)*6*npart);
        yprime(ip) = Ctot(ip+5*npart+(islice-1)*6*npart);
        Ctot(ip+npart+(islice-1)*6*npart) = phiabs_fin(ip);
        Ctot(ip+2*npart+(islice-1)*6*npart) = xpos(ip);
    end
end

bunchpartfile = strcat(dirname,'mg.out.dpa');
Cpoint = fopen(bunchpartfile,'w+');
Ccount = fwrite(Cpoint,Ctot,'double');
fclose(Cpoint);
%%
figure(31)
plot(phiabs,gamma,'b+')
hold on
plot(phiabs_fin,gamma,'go')
hold off
bcomplex0=mean(mean(exp(1i*phiabs)));
bcomplex=mean(mean(exp(1i*phiabs_fin)));
titlestr=sprintf('R56=%.2e phaseshift=%.2f bfactor_i=%.2f bfactor_f=%.2f',R56buncher,phaseshift,abs(bcomplex0),abs(bcomplex));
title(titlestr);
figure(32)
plot(xpos,gamma,'go');


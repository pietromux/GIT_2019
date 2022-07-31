%plot_lps_td plots the entire longitudinal phase space of the beam 

Cpoint = fopen(outputpartfile,'r');
Ctot = fread(Cpoint,6*npart*niter,'double');
fclose(Cpoint);
for islice = 1:niter
    for ip =1:npart
        gamma_sp(islice,ip) = Ctot(ip+(islice-1)*6*npart); 
        phiabs_sp(islice,ip) = Ctot(ip+npart+(islice-1)*6*npart)+islice*2*pi;
        xpos(islice,ip) = Ctot(ip+2*npart+(islice-1)*6*npart);
        xprime(islice,ip) = Ctot(ip+4*npart+(islice-1)*6*npart);
        ypos(islice,ip) = Ctot(ip+3*npart+(islice-1)*6*npart);
        yprime(islice,ip) = Ctot(ip+5*npart+(islice-1)*6*npart);
        phiplot(islice,ip) = mod(phiabs_sp(islice,ip),2*pi);
    end
end
figure(1)
plot(phiabs_sp,gamma_sp,'.')
for islice=curpeakpos
     gamma_p=gamma_sp(islice,:);
     phiplot_p = phiplot(islice,:);
     theta = -pi:0.001:pi;
     psi0=pi/4; 
     theta2=pi-psi0;
     gammares = sqrt(git_par.lambda_w0*(1+aw0^2)/git_par.lambda/2);
     bucket=sqrt(git_par.lambda_w0/git_par.lambda*K_laser(end)*aw0*(cos(theta)+theta*sin(psi0)-cos(theta2)-theta2*sin(psi0)));
     bucketplot=-bucket( bucket>0)+gammares; %y-value of lower bucket
     bucketplot2=bucket( bucket>0)+gammares; %y-value of upper bucket
     %thetaplot=theta( bucket>0)+median(phiplot)-std(phiplot)*0.7+phaseadd(phasej); %domain of bucket
     figure(2)
     scatter(phiplot_p,gamma_p,'.')   %%plot the particle energy vs phas
     hold on
     plot(theta(bucket>0),bucketplot,'.k') %plot lower bucket
     plot(theta(bucket>0),bucketplot2,'.k') %plot upper bucket
     hold off     
end

% bcomplex=mean(mean(exp(1i*phiabs_sp)));
% bfactor=abs(bcomplex);
% bphase=angle(bcomplex);

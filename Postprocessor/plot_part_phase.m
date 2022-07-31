%plot_part_phase is specified to make 6 subplots of particle energy vs. phase for time-indep GIT run
%of TESSA266 that output around 35GW power. 

if exist('num','var')==0
    num=77;
    figure(num)
set(gcf,'name','GIT phase plot');
end
if exist('phasej','var')==0
    phasej=1;
    num=77;
    
end
if phasej>6
    phasej=1;
    num=num+1;
end
GITphasefig=figure(num);
    
% subplot(2,3,phasej);

%%Make the bucket for TESSA266; all parameters are specified for t-indep
%%run that outputs around 35GW. 
%%phasej==1 is after running prebuncher and before R56
%%phasej==2 is after running prebuncher R56 (prop_part)
%%phasej==3 to 6 are after each undulator sections. 

%some change
Cpoint = fopen(outputpartfile,'r');
Ctot = fread(Cpoint,6*npart*niter,'double');
fclose(Cpoint);
for islice = 1:niter
    for ip =1:npart
        gamma_sp(islice,ip) = Ctot(ip+(islice-1)*6*npart); 
        phiabs_sp(islice,ip) = Ctot(ip+npart+(islice-1)*6*npart);
        xpos(islice,ip) = Ctot(ip+2*npart+(islice-1)*6*npart);
        xprime(islice,ip) = Ctot(ip+4*npart+(islice-1)*6*npart);
        ypos(islice,ip) = Ctot(ip+3*npart+(islice-1)*6*npart);
        yprime(islice,ip) = Ctot(ip+5*npart+(islice-1)*6*npart);
    end
end
for islice=1:niter
    phiabs=phiabs_sp(islice,:);
    gamma=gamma_sp(islice,:);
phiplot = mod(phiabs,2*pi);

theta = -pi:0.001:pi;
psi0=pi/4-0.0; %%define a phase that produce reasonable bucket
theta2=pi-psi0;
gamlim=[0,737,737,747,750,750]; %%define energy limit that produce reasonable bucket
phaseadd=[0,0,0,0,0,-pi/7]; %%some plots are shifted..

bucket=sqrt(git_par.lambda_w0/git_par.lambda*K_laser(end)*aw0*(cos(theta)+theta*sin(psi0)-cos(theta2)-theta2*sin(psi0)));

buck_phiabs = sqrt(git_par.lambda_w0/git_par.lambda*K_laser(end)*aw0*(cos(phiabs)+mod(phiabs,2*pi)*sin(psi0)-cos(theta2)-theta2*sin(psi0)));
buck_gam = gamma-mean(gamma);
sum(-buck_phiabs<buck_gam & buck_gam<buck_phiabs & phiplot < 2.35 & phiplot > -0.2)
figure(35)
plot(phiplot,buck_phiabs,'r.')
hold on
plot(phiplot,-buck_phiabs,'g.')
plot(phiplot,buck_gam,'bo')
hold off
figure
bucketplot=-bucket( bucket>0)+median(gamma(gamma<gamlim(phasej))); %y-value of lower bucket
bucketplot2=bucket( bucket>0)+median(gamma(gamma<gamlim(phasej))); %y-value of upper bucket
thetaplot=theta( bucket>0)+median(phiplot)-std(phiplot)*0.7+phaseadd(phasej); %domain of bucket
if islice==round(niter/2) %only plot for one of the slices
    scatter(phiplot,gamma,40,cool(length(gamma))-.1*islice/niter,'.')   %%plot the particle energy vs phas
    hold on
plot(thetaplot,bucketplot,'.k') %plot lower bucket
plot(thetaplot,bucketplot2,'.k') %plot upper bucket
end
end
hold off

bcomplex=mean(mean(exp(1i*phiabs_sp)));
bfactor=abs(bcomplex);
bphase=angle(bcomplex);


phasej=phasej+1; %increase phasej for each run.
% 
bstring=sprintf('bunch=%.3f z=%.2f',bfactor,z(j)); %outputs the lastest bunch factor and z
title(bstring,'FontSize',8);
xlim([0 2*pi]);
xlabel('\phi')
ylabel('\gamma')
set(gca,'xtick',0:pi/2:2*pi,'xticklabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
% title(sprintf('phase space at z=%.3f for %s',z(j),mat2str(sep)));
% figname=sprintf('phase z=%.4f sep=%s',z(j),mat2str(sep0));


%One could save the phases by these settings
% figfiledir=[datadir,'phases\',datetimenow,'\'];
% if exist(figfiledir,'dir') == 0
%     mkdir(figfiledir)
% end
% figfilename=[figfiledir,figname,'.png'];
% saveas(gcf,figfilename); 
% set(gcf, 'PaperUnits', 'inches');set(gcf, 'Position', [1 1 1500 1200]); set(gcf, 'PaperPosition', [0 0 12.5 6.5]);

outfig=figure(501);
clf(gcf);
set(gcf,'Name','slideplot2_power')
set(gca, 'Color', 'none')

subplot(2,3,1)
plot(z,power,'.');  
ylabel('power(W)');
xlabel('z(m)');
xlim([0 max(z)])
ylim([0 max(power)])
hold on

%enhance_plot
title(sprintf(' max=%.3e W at z= %.3f ',power(end),z(end)),'FontSize',15)

subplot(2,3,3)
set(gca, 'Color', 'none')

hold on
plot(z(2:end),xrms(2:end),'r','DisplayName','xrms')
plot(z(2:end),yrms(2:end),'k','DisplayName','yrms')
legend('show')
hold off
avgsigx=mean(xrms(2:end));
avgsigy=mean(yrms(2:end));
avgsigxy=sqrt(avgsigx*avgsigy);
ylabel('beam size (m)')
xlabel('z(m)')
title(sprintf('<\\sigma>=%.2e',avgsigxy))
xlim([0 max(z)])
ylim([min([xrms,yrms]) max([xrms,yrms])]);

%enhance_plot

subplot(2,3,2)
set(gca, 'Color', 'none')


plot(z(2:end),energy(2:end));
title('Beam energy')
ylabel('gamma')
xlabel('z(m)')
%enhance_plot
set(gca, 'Color', 'none')
xlim([0 max(z)])

subplot(2,3,4)
set(gca, 'Color', 'none')

plot(z,K_undulator(1:length(z)),'.')
title('K undulator')
ylabel('K und')
xlabel('z(m)')
%enhance_plot
xlim([0 max(z)])
% ylim([min(K_undulator(K_undulator~=0)) max(K_undulator)])


subplot(2,3,5)
set(gca, 'Color', 'none')

plot(z,bunching)
xlim([0 max(z)])

ylabel('bunching')
xlabel('z(m)')
hold on
title('bunching')
%enhance_plot

subplot(2,3,6)
plot(z,r_size,'g','DisplayName','rsize')
hold on
plot(z,xrms,'r','DisplayName','xrms')
plot(z,yrms,'k','DisplayName','yrms')
hold off
legend('show')
ylabel('spot size (m)')
xlabel('z(m)')
xlim([0 max(z)])

%enhance_plot
title('Beam and radiation rms spot size','FontSize',15)
set(gca, 'Color', 'none')
figure(5)
subplot(4,3,6)
plot(z,p_mid)
title('On axis intensity')
ylabel('Intensity')
xlabel('z(m)')

% subplot(4,3,8)
% hold on
% plot(z,xrms,'r','DisplayName','xrms')
% plot(z,yrms,'k','DisplayName','yrms')
% % avgsigx=trapz(z,xrms)/(z(length(z))-z(1));
% % avgsigy=trapz(z,yrms)/(z(length(z))-z(1));
% ylabel('beam size (m)')
% xlabel('z(m)')
%     
% 
% hold off
% legend('show')
% 
% title('Beam size')
% 
% outres.z= z;
% outres.K_undulator = K_undulator;
% outres.lambda_undulator = lambda_undulator;
% outres.power = power;
% outres.phislice = phislice;
% outres.energy = energyslice;
% outres.bunching = bunching;
% outres.trapped = ncslice;
% outres.p_mid = p_mid;
% outres.r_size = r_size;
% 
% outres.xrms = xrms;
% outres.d3fac = d3fac;
% outres.resphase = resphase;
% outres.gammaresonant = gammaresonant;
% outres.theta = thetaslice;
% outres.deltatheta = deltatheta;
% 
% subplot(4,3,9)
% plot(z,bunchphases)
% title('bunchphase vs. z')
% %if(~itdp)
% %outres.zwg = zwg;
% %outres.powerwg = powerwg;
% %outres.bunchingwg = bunchingwg;
% %outres.xrmswg = xrmswg;
% %outres.r_sizewg = r_sizewg;
% %outres.phiwg = phiwg;
% %end
% 
% %end
% avgsigtxt=[sprintf(' lambdau= %.5f\n aw0 = %.3f\n x = %s \n sep = %s \n sep_w = %s \n B0 = %.2f T/m \n z(end) = %.3f \n power(end) =%.3e\n rxbeam = %.2e\n rybeam = %.2e\n alphax = %.3f\n alphay = %.3f\n avgsigx = %.2e\n avgsigy = %.2e\n avgsigxy = %.2e\n zrayl=%.2f\n zwaist=%.2f\n bunchphase=%.2f\n deltapsi=%.2f\n psi_0=%.2f\n',lambda_w0, git_par.aw00,mat2str(x),mat2str(sep),mat2str(sep/lambda_w0),git_par.quadgradient,z(end),power(end),rxbeam,rybeam,alphax,alphay,avgsigx,avgsigy,sqrt(avgsigx*avgsigy),git_par.zrayl,git_par.zwaist,git_par.bunchphase,git_par.deltapsi,git_par.psi_0)];
% % avgsigtxt=[sprintf('bunch =%.2f \n bunchphase =%.2f \n power =%.3e \n' , git_par.bunch,git_par.bunchphase,power(end))]
% annotation('textbox',[0.01,.7,.1,.1],'string',avgsigtxt);
% dir='Y:\youna\matlab\git\gitfig\32\';
set(gcf, 'PaperUnits', 'inches');
set(gcf,'Units','inches');
set(gcf, 'Position', [1 1 13 6])
set(gcf, 'PaperPosition', [0 0 13 6]);

for nn=1:1:6
   
    subplot(2,3,nn)
    set(gca, 'Color', 'none')
end
    pmax0=power(end);
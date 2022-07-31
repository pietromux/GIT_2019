zpos=slicenum*git_par.lambda; timearray=slicenum*xlamds/3e8;

outfig=figure(501);
clf(gcf);
subplot(4,3,1);
set(gcf,'name','tdep power plot')
hold on
powerarray=power;
Msize=size(power);
zarray=z;
pmaxend=max(max(power));
ColorSet=colormap(cool(Msize(2)));
for jp=1:Msize(2)
%     ColorSet(jp,:)=[1-slicenum(jp)/slicenum(end) 0 slicenum(jp)/slicenum(end)];
plot(slicenum,powerarray(:,jp)/pmaxend,'color',ColorSet(jp,:))
end
plot(slicenum,current/max(current),'k.-');
hold on
plot(slicenum,p_rad_0/max(p_rad_0),'r.-');
plot(slicenum,p_rad_1/max(p_rad_1),'g.-');
plot(slicenum,power(:,2)./max(power(:,2)),'b.-');
set(gcf, 'Colormap', ColorSet);
caxis([z(1) z(end)]);
xlabel('slicenum');
ylabel('power(/max) and current (/max)');
cb=colorbar;
ylabel(cb,'distance in undulator z(m)')
title('k=current, r=seed, g=after pb before slip, b=und 1st period')

subplot(4,3,2);
psize=size(power);
colorset=cool(psize(1));
for jp=1:psize(1)
    hold on
plot(z,power(jp,:)/1e9,'Color',colorset(jp,:))
end
xlabel('z(m)')
ylabel('power(GW)')
set(gcf,'Colormap',colorset);
cb=colorbar;
ylabel(cb,'time(ps)')
ylim(cb,[min(timearray*1e12) max(timearray*1e12)])
subplot(4,3,3);
[C,h]=contourf(slicenum*266e-9/3e8/1e-12,z,power'/pmaxend);
hold on
plot(slicenum*266e-9/3e8/1e-12,current/max(current),'k.-')
set(h,'LineColor','none')
colormap('cool')
cb=colorbar;  
ylabel(cb,'power(/max)')
xlabel('time(ps)')
ylabel('z(m)')

subplot(4,3,6)
colorsize=size(energy);
colorset=cool(colorsize(1));
for jp=1:colorsize(1)
    hold on
plot(z,energy(jp,:),'Color',colorset(jp,:))
end
set(gcf,'Colormap',colorset);
cb=colorbar;
ylabel(cb,'time(ps)')
ylim(cb,[min(timearray*1e12) max(timearray*1e12)])
ylabel('\gamma')
xlabel('z(m)')
title('Beam energy')
subplot(4,3,9)
colorsize=size(bunching);
colorset=cool(colorsize(1));
for jp=1:colorsize(1)
    hold on
plot(z,bunching(jp,:),'Color',colorset(jp,:))
end
meanbunching=0;
for jp=1:size(power,2)
    meanbunching(jp)=sum(bunching(:,jp).*power(:,jp))/sum(power(:,jp));
end
plot(z,meanbunching,'k-','LineWidth',2);
xlabel('z(m)')
ylabel('bunching')
set(gcf,'Colormap',colorset);
cb=colorbar;
ylabel(cb,'time(ps)')
ylim(cb,[min(timearray*1e12) max(timearray*1e12)])
title('bunching')
subplot(4,3,7)
hold on
plot(z,p_mid)
ylabel('|E|^2')
xlabel('z(m)')
title('On axis intensity')
cb=colorbar;
ylabel(cb,'time(ps)')
ylim(cb,[min(timearray*1e12) max(timearray*1e12)])
subplot(4,3,8)
colorsize=size(r_size);
meanrsize=[];
for jp=1:colorsize(1)
    hold on
plot(z,r_size(jp,:),'Color',colorset(jp,:))
end
for jp=1:size(power,2)
meanrsize(jp)=sum(r_size(:,jp).*power(:,jp))/sum(power(:,jp));
end
plot(z,meanrsize,'k-','LineWidth',2)
set(gcf,'Colormap',colorset);
cb=colorbar;
ylabel(cb,'time(ps)')
ylim(cb,[min(timearray*1e12) max(timearray*1e12)])
hold on
plot(z,xrms,'b','DisplayName','xrms')
plot(z,yrms,'r','DisplayName','yrms')
hold off
rsizetitle=sprintf('b=xrms, r=yrms, color=rsize, rsize min=%.2f um',min(min(r_size))*1e6); 
title(rsizetitle);
xlabel('z(m)')
ylabel('spot size(m)')

powermaxend=max(power(:,end));
powerinitial=375e6*max(current);
powereff=powermaxend/powerinitial*100;

powerend=power(:,end);
zpos2=lambda*(1:1:size(power,1));
Uinitial=trapz(zpos./3e8,current)*gamma0*0.511e6;
Ufinal=trapz(zpos2./3e8,powerend);
Ueff=Ufinal/Uinitial*100;

for j=1:size(power,2)
powermaxz(j)=max(power(:,j));
powereffz(j)=powermaxz(j)/powerinitial*100;
powerz_weighted(j)=sum(power(:,j).*power(:,j))/sum(power(:,j));
ueffz(j)=trapz(zpos./3e8,power(:,j))/Uinitial*100;
ueffz2(j)=sum(power(power(:,j)>9e9,j))*xlamds/3e8/Uinitial*100;
ufinalz(j)=sum(power(power(:,j)>9e9,j))*xlamds/3e8;
end
subplot(4,3,10)
plot(z,K_undulator)
xlabel('z(m)')
ylabel('K')
% ylim([min(K_undulator(K_undulator~=0)) max(K_undulator)]);
title('tapering')

subplot(4,3,11)
plot(z,xrms,'b')
hold on
plot(z,yrms,'r')
xlabel('z(m)')
ylabel('e-beam spots size')
title('e-beam spot size, b=xrms r=yrms')
%plot_spectral_data_contour_sub;


%streff=sprintf(' powermax= %.2E W\n Ueff_1 = %.2f %%, Ueff_5 = %.2f %%\n Ueff_{all} = %.2f %%\n powereff = %.2f %% \n maxcurrent= %.2f kA \n curlen(RMS)= %.2f um, %.f slices',powermaxend,ueff_spec_1,ueff_spec_5,Ueff,powereff,max(current)/1e3,git_par.curlen*1e6,git_par.curlen/git_par.lambda);
%annotation('textbox',[.42,.83,.1,.1],'string',streff,'FontWeight','bold');

% selectj=[2,5,6,8,9];
% for jjs=1:length(selectj)
%     subplot(3,3,selectj(jjs));
%     xlim([0 z(end)+.2])
% end

   
set(gcf, 'PaperUnits', 'inches');
set(gcf,'Units','inches');
set(gcf, 'Position', [1 1 20 8])
set(gcf, 'PaperPosition', [0 0 20 8]);


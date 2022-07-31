% Plot spectral data
xlamd=git_par.lambda_w0;
xlamds=git_par.lambda;
clear zlocations pulseenergy
outfig=figure(501);
c=3e8;q_e=1.6e-19;eps0=4*pi*1e-07;
colorset=cool(psize(2));
kw=2*pi/xlamd;
m_e=9.11e-31;
hbar=6.582e-16;
fundpower=[];sidebandpower=[];
zlocations=z;
zindices=ones(1,length(z));
tvector=slicenum*xlamds/c;
specwind=10;
p_spec=[];
p_spec_all=[];
  for n=1:length(zlocations)
      legendstring=['z = ',num2str(zlocations(n)),'m'];
      poweramp=power(:,n);
      ephase=phi_mid_M(:,n);
      powerfield=sqrt(poweramp).*exp(1i*ephase);
      [specvar,omegavar]=g3spectrum2(powerfield,xlamds,zsep);
 specfiltered=specvar;specfiltered2=specvar;
 omegamax=1e-3; omegamin=-1e-3;
 sidebandindex=omegavar>omegamax & omegavar < omegamin;
 specfiltered(sidebandindex)=0;
 fundpower(n)=trapz(specfiltered)/trapz(specvar);
 sidebandpower(n)=(trapz(specvar)-trapz(specfiltered))/trapz(specvar);
    pulseenergy(n)=trapz(tvector,power(:,n));
subplot(4,3,4)
   hold on
farray=((1:1:nslice)-nslice/2)/nslice*c/xlamds+c/xlamds;
wavelenvar=c./farray;
[minval,minind]=min(abs(wavelenvar-xlamds));
for jw=1:specwind
p_spec(n,jw)=sum(specvar(minind-jw:minind+jw))*xlamds/3e8/length(specvar);
end
plot(wavelenvar,specvar,'-','Color',colorset(n,:));

p_spec_all(n)=sum(specvar)*xlamds/3e8/length(specvar);
  end
dlam=wavelenvar(minind-1)-wavelenvar(minind);
title('energy spectral density')
     xlabel('Wavelength(m)')
     ylabel('P (\omega) [arb. units]','FontSize',16)
     cb=colorbar;
     colormap(cb,cool);
     ylabel(cb,'z(m)')
     caxis([zlocations(1) zlocations(end)]);
          xlim([266e-9-5e-9 266e-9+5e-9])

% legend(legendstring); legend boxoff
subplot(4,3,5)
colorlist=cool(specwind);
for jw=1:specwind
    plot(zlocations,p_spec(:,jw),'Color',colorlist(jw,:));
    hold on
end
    plot(zlocations,p_spec_all,'-k');

    xlabel('z(m)')
    ylabel('energy(J)')
    wavelenstr=sprintf('ueff near \\lambda_s \\pm 1-10 d\\lambda(color) and all (black) where d\\lambda=%.1fnm',dlam*1e9);
title(wavelenstr);
ueff_spec=p_spec/Uinitial;
ueff_spec_all=p_spec_all/Uinitial;
ueff_spec_5=ueff_spec(end,5)*100;
ueff_spec_1=ueff_spec(end,1)*100;
cb=colorbar;
caxis([1 10])
ylabel(cb,'\lambda \pm n d\lambda');

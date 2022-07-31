% Plot spectral data
function p_spec_all=plot_spec(powerfield,git_par)
xlamd=git_par.lambda_w0;
xlamds=git_par.lambda;
zsep=git_par.zsep;
nslice=git_par.nslice;
c=3e8;q_e=1.6e-19;eps0=4*pi*1e-07;
kw=2*pi/xlamd;
m_e=9.1e-31;
hbar=6.582e-16;
fundpower=[];sidebandpower=[];
slicenum=1:size(powerfield,1);
tvector=slicenum*xlamds/c;
specwind=10;    
p_spec=[];
p_spec_all=[];
[specvar,omegavar]=g3spectrum2(powerfield,xlamds,zsep);
specfiltered=specvar;specfiltered2=specvar;
omegamax=1e-3; omegamin=-1e-3;
sidebandindex=omegavar>omegamax & omegavar < omegamin;
specfiltered(sidebandindex)=0;
hold on
farray=((1:1:nslice)-nslice/2)/nslice*c/xlamds+c/xlamds;
wavelenvar=c./farray;
[minval,minind]=min(abs(wavelenvar-xlamds));
for jw=1:specwind
p_spec(jw)=sum(specvar(minind-jw:minind+jw))*xlamds/3e8/length(specvar);
end
plot(wavelenvar,specvar)
p_spec_all=sum(specvar)*xlamds/3e8/length(specvar);
title('energy spectral density')
xlabel('Wavelength(m)')
ylabel('P (\omega) [arb. units]','FontSize',16)
xlim([266e-9-5e-9 266e-9+5e-9])
%function [ propfield ] = free_prop_spectr( git_par, M, Mag, field )

%Spectral fourier Wavefront Propagation algorithm
% Inputs 
% field  = Complex input field in space domain
% dx space step of the field array
% M propagation matrix
% Mag magnification of final grid
%%%%%%%
% Output
% propfield = complex output field in space domain
% only works for single slice

%% Load field
Ltot=15;

outputfieldfile = strcat(git_par.dirname,'mg.out.dfl');
Cfieldfile = fopen(outputfieldfile);
Cfield = fread(Cfieldfile,2*git_par.ncar*git_par.ncar*niter,'double');
fclose(Cfieldfile);                                     

zw=1.4;
Lu=6.304;
f=git_par.f;
d=git_par.d;
L1=Ltot/2-(Lu-zw)+d;
L2=Ltot/2-zw-d;
poweramp=[];poweramp_out=[];rms_in=[];rms_out=[];
for islice = 1:niter
    iindex = git_par.ncar*git_par.ncar*(islice-1)*2;
    recfield = Cfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1);
    imcfield = Cfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar);

    dx = 2*git_par.dgrid/git_par.ncar;
    lambda = git_par.lambda;
    k = 2*pi/lambda;
    field = reshape((recfield-1i*imcfield),ncar,ncar); %note that by convention we need the complex coniugate of Genesis field
  
atten = 1;


lens = [1 0 ; -1/f 1]; 
drift1 = [1 L1; 0 1];
drift2= [1 L2; 0 1];
M = drift2*lens*drift1;
Mag =1;

    A = M(1,1);
    B = M(1,2);
    C = M(2,1);
    D = M(2,2);

    [Mx, My] = size(field);
    a1 = dx*(Mx+1)/2;
    a3 = a1;
    a2 = Mag*a1;
    a4 = a2;

    dkx = 2*pi/(k*dx*Mx) ; dky = 2*pi/(k*dx*My); 
    nx = ((1:Mx)-(Mx+1)/2);
    ny = ((1:My)-(My+1)/2);
    [kx,ky] = meshgrid(nx*dkx,ny*dky);
    [x1,y1] = meshgrid(nx*dx,ny*dx);
    [x2,y2] = meshgrid(nx*dx*Mag,ny*dx*Mag);

    vfield = field.*exp(-1i*(pi*(A-Mag).*(x1.^2+y1.^2)/B/lambda));
    Fresnel = Mag/B/lambda;
    psift = (fftshift(fft2(vfield))); 
    psi_k = psift.*exp(1i.*(kx.^2+ky.^2)*k^2/4/pi/Fresnel);

    propvfield = (ifft2(ifftshift(psi_k)));
    propfield = propvfield.*exp(-1i*(pi*(D-1/Mag).*(x2.^2+y2.^2)/B/lambda));
    phase_out = angle(propfield((ncar+1)/2,(ncar+1)/2));
    propfield = propfield.*exp(-1i*phase_out).*atten;
    %% output
    power_in = sum(sum(abs(field).^2));
    power_out= sum(sum(abs(propfield).^2));
    rms_in(islice) = 1.414*sqrt(sum(sum(abs(field).^2.*(x1.^2+y1.^2)))/power_in);
    rms_out(islice) = 1.414*sqrt(sum(sum(abs(propfield).^2.*(x2.^2+y2.^2)))/power_out);
    phase_out = angle(propfield((ncar+1)/2,(ncar+1)/2));
    phase_in=angle(field((ncar+1)/2,(ncar+1)/2));

% formatSpec = 'Rad_propagation. Power %.3e ->  %.3e.   \n rms %.3e -> %.3e . Phase %.3e \n';
% fprintf(formatSpec,power_in,power_out,rms_in,rms_out,phase_out);
% end
    repropfield = real(propfield);
    impropfield = imag(propfield);
    outrefield =reshape(repropfield,ncar*ncar,1);
    outimfield =reshape(impropfield,ncar*ncar,1);
    outfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1) = outrefield;
    outfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar) = -outimfield;
    poweramp(islice)=power_in;
    poweramp_out(islice)=power_out;
end
    
titlestr=sprintf('L=%.2f F=%.2e rmsin=%.2e m rmsout=%.2e m',L1,f,mean(rms_in(isnan(rms_in)==0)),mean(rms_out(isnan(rms_out)==0)));

%output the contour plot of the maximum amplitude slice
[maxval,maxind]=max(poweramp);
iindex = git_par.ncar*git_par.ncar*(maxind-1)*2;
refield_m=Cfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1);
imfield_m=Cfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar);
field_m=refield_m-1i*imfield_m;
outrefield_m=outfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1);
outimfield_m=outfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar);
propfield_m=outrefield_m-1i*outimfield_m;

recircfig1=figure(22);
clf
set(gcf,'name','recirc check on power and rms')
subplot(2,1,1)
plot(poweramp,'*')
hold on
plot(poweramp_out)
legend('poweramp in','poweramp out')
subplot(2,1,2)
plot(rms_in)
hold on
plot(rms_out)
legend('rmsin','rmsout')

recircfig2=figure(24);
clf;
set(gcf,'name','prop_spectr_recirc contour')
subplot(2,1,1)
contourf(abs(reshape(field_m,ncar,ncar)))
title(titlestr);
legend('input field')
subplot(2,1,2)
contourf(abs(reshape(propfield_m,ncar,ncar)))
legend('output field')

%saveas(gcf,[figdir,'propspec_contour_',num2str(iss),'.png']);

recircfig3=figure(25);
clf;
set(gcf,'name','intput and output field')
plot(rms_in)
hold on
title('input field')
plot(rms_out)
hold off
title('output field')
legend('rmsin','rmsout')
title(titlestr)

if saveplots==1
saveas(recircfig1,[figdir,'propspec_power_',num2str(iss),'.png']);
saveas(recircfig2,[figdir,'propspec_contour_',num2str(iss),'.png']);
saveas(recircfig3,[figdir,'propspec_rms_',num2str(iss),'.png']);
end
% 
% 
% 
outputfieldfile = strcat(git_par.dirname,'mg.out.dfl');

Cfieldfile = fopen(outputfieldfile,'w+');
Ccount = fwrite(Cfieldfile,outfield','double');
fclose(Cfieldfile);

% % 
% % % 
% 

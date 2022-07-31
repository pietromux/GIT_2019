% TESSA amplifier Genesis Informed Tapering Scheme 
% P. Musumeci 2019 version 

%%%%------ This script should never be modified. --------------%%%%
%%%%------ In order to change parameters , change Git_par.m file only ------%%%%
%%%%------ It runs Genesis for NWIG periods multiple times. Undulator
%%%%             period and K can vary from period to period. 

%% Initialize local variables
helical = git_par.helical;             
period_tapering = git_par.period_tapering; 
verbose = git_par.verbose;             
outputfrequency = git_par.outputfrequency; 
correct_theta = git_par.correct_theta;                 % Use for phase-shifters
correct_R56 = git_par.correct_R56;

fac0 = git_par.fac0;               
psi_0 = git_par.psi_0;
quadrupole_lattice = git_par.quadrupole_lattice; 
dircyg = git_par.dircyg;
dirname=git_par.dirname;
homedir=git_par.homedir;

formatSpec = 'TESSA simulation. fac0 %.3f psi_0 %.3f \n';
fprintf(formatSpec,fac0,psi_0);

z_stop = git_par.z_stop;          % Simulation length

itgaus = git_par.itgaus;     % 0 parabolic ; 1 gaussian ; 2 uniform transverse profiles
iscan = git_par.iscan;
nscan = git_par.nscan;
svar =  git_par.svar;
delz = git_par.delz;
itdp = git_par.itdp;
iotail = git_par.iotail;
nslice = git_par.nslice;
zsep = git_par.zsep;

slippage = git_par.slippage;
slippagetot = git_par.slippagetot;

if (itdp)
    niter = nslice;
else
    niter = nscan;
end
wcoefz = git_par.wcoefz;
lambda_w0 = git_par.lambda_w0;
und_type = helical;
npart = git_par.npart;
gamma0 = git_par.gamma0; 
delgamma = git_par.delgamma;
rxbeam = git_par.rxbeam;
rybeam = git_par.rybeam;
emit = git_par.emit;
curpeak = git_par.curpeak;
bunch = git_par.bunch;
bunchphase = git_par.bunchphase; 
prad0 = git_par.prad0;
lambda = git_par.lambda;
zrayl = git_par.zrayl;
zwaist = git_par.zwaist;
nwig = git_par.nwig;
ncar = git_par.ncar;
dgrid = git_par.dgrid;
w0 = sqrt(lambda*zrayl/pi);
went = w0*sqrt(1+zwaist^2/zrayl^2);
if correct_theta 
deltatheta=git_par.deltatheta;
end

%% Define period and undulator K in 3 cases.  
% period tapering > 0. Uses Halbach undulator builder equation
% period tapering = 0. Keeps constant period
% period tapering < 0. Pre-assigned values 

if(period_tapering>0)
    gap = git_par.gap;
    Brem = git_par.Brem;
    B_und = @(eta) 1.6*Brem*exp(-pi*gap/eta);
    K_und = @(eta) B_und(eta)*eta*93;
    dK_und =@(eta) 93*B_und(eta)+93*eta*(1.6*pi*Brem*gap)/(eta^2*exp((pi*gap)/eta));
    if(~helical)
        res =@(eta) eta*(1+K_und(eta)^2/2)/2/gamma0.^2-lambda;
        lambda_w = fzero(res,lambda_w0);
        aw0 = K_und(lambda_w)/1.414;
    else
        rest =@(eta) eta*(1+K_und(eta)^2)/2/gamma0.^2-lambda;
        lambda_w = fzero(rest,lambda_w0);
        aw0 = K_und(lambda_w);
    end
elseif (period_tapering ==0)
    lambda_w = lambda_w0;
    aw0 = sqrt(lambda*2*gamma0.^2/lambda_w - 1); 
elseif (period_tapering==-1)
    lambda_w = Nocibur_periods(1);
    aw0 = Nocibur_K_init(1);
end

k_w = 2*pi/lambda_w;

%% Use constant beam size (quadrupole_lattice = 0) or simulate real quadrupole lattice (quadrupole_lattice = 1). 
% Note that if drift regions are to be simulated, one must use a real quadrupole lattice
% (i.e. quadrupole_lattice = 1) since beta-function is not constant

if quadrupole_lattice
    rxbeam = git_par.rxbeam_match;
    rybeam = git_par.rybeam_match;
    alphax = git_par.alphax;
    alphay = git_par.alphay;
    xkx = git_par.xkx;
    xky = git_par.xky;
else
     beta_function = rxbeam^2*gamma0/emit;
     xkx = emit^2/(rxbeam^4*aw0^2*k_w^2);
     xky = xkx;
     alphax = 0;
     alphay = 0;
     drift = 0;
end

%% Initialize file names and output variables
filename = strcat(dirname,'mg.in');
inputfilename = 'TESS3.0/mg.in';
outputfilename = 'TESS3.0/mg.out';
radfilename = 'TESS3.0/mg.rad.in';
beamfilename = 'TESS3.0/mg.beam.in';

clear z Epart0 K_undulator
j = 1;
incK = 0;
incl = 0;
z(1) = 0;

if (~git_par.beamfile & ~prebunched)
    bunching(1:niter,1) = bunch;
    energy(1:niter,1) = gamma0;
    xrms(1:niter,1) = rxbeam;
    yrms(1:niter,1) = rybeam;
elseif (prebunched)
    bunching(1:niter,1) = bunching(1:niter,end); bunching(:,2:end) = [];
    xrms(1:niter,1) = xrms(1:niter,end);    xrms(:,2:end) = [];
    yrms(1:niter,1) = yrms(1:niter,end);    yrms(:,2:end) = [];
    energy(1:niter,1) = energy(1:niter,end);     energy(:,2:end) = [];    
    energy
end    

if (~git_par.radfile & ~prebunched)
    r_size(1:niter,1) = went/2;
    power(1:niter,1) = prad0;
elseif(prebunched)
    r_size(1:niter,1) = r_size(1:niter,end); r_size(:,2:end) = [];
    power(1:niter,1) = power(1:niter,end); power(:,2:end) = [];
    power
end
lambda_undulator(1) = lambda_w;
K_undulator(1) = aw0;
xkx0 = xkx;
offset = git_par.offset;
awd=git_par.awd;
%% Main simulation loop 
while (z(end) < z_stop)
    tstart = tic;
%% Magnetic lattice implementation.     
    if(quadrupole_lattice==1)
        drift = 1;
        if (lat(1,j))
            lat(1,j) = aw0;
            correct_R56 = 0;
            offset = git_par.offset;
            iotail = git_par.iotail;
            awd=0;
        else
            lat(1,j) = 0;
%         lat(3,j) = aw0;
            lat(3,j) = git_par.awd;
            correct_R56 = git_par.correct_R56;
            offset = 0;
            iotail = 1;
%         awd=lat(3,j);
            awd=git_par.awd;
        end
        
        if j>2 & lat(1,j)~=0 & lat(1,j-1)==0
            at_undentrance=1;
            at_undexit=0;
        elseif j+1<z_stop/lambda_w & lat(1,j)~=0 & lat(1,j+1)==0
            at_undexit=1;
            at_undentrance=0;            
        else
            at_undentrance=0;
            at_undexit=0;
        end
        
        magfileName=strcat(dirname,'mg.mag.in');
        fid=fopen(magfileName,'w');   
        fprintf(fid,'%s\n','? VERSION=1.0');
        fprintf(fid,'%s %f\n','? UNITLENGTH=',lambda_w);
        fprintf(fid,'%s %f %.f %.f\n', 'AW', lat(1,j), 1, 0);
        fprintf(fid,'%s %f %.f %.f\n', 'QF', lat(2,j), 1, 0);
%     fprintf(fid,'%s %f %.f %.f\n', 'AD', lat(3,j), 1, 0);
        fprintf(fid,'%s %f %.f %.f\n', 'AD', awd, 1, 0);
        fclose(fid);
     elseif quadrupole_lattice==2
        correct_R56=git_par.correct_R56;
        drift = git_par.drift;
    end
    
    %% Write input file for Genesis
fid = fopen(filename,'wt');
fprintf(fid, ' $newrun \n');
fprintf(fid, ' aw0 = %E \n', aw0);
fprintf(fid, ' xkx = %E \n', xkx);
fprintf(fid, ' xky = %E \n', xky);
fprintf(fid, ' wcoefz = %E %E %E\n', wcoefz(1), wcoefz(2),wcoefz(3));
fprintf(fid, ' xlamd = %E \n', lambda_w);
fprintf(fid, ' iwityp = %d \n', und_type);
fprintf(fid, ' awd = %E \n', 1); % AD in magin file overrides awd
fprintf(fid, ' npart = %d \n', npart);
fprintf(fid, ' gamma0 = %E\n', gamma0);
fprintf(fid, ' delgam = %E\n', delgamma);
fprintf(fid, ' rxbeam = %E\n', rxbeam);
fprintf(fid, ' rybeam = %E\n', rybeam);
fprintf(fid, ' alphax = %E\n', alphax);
fprintf(fid, ' alphay = %E\n', alphay);
fprintf(fid, ' emitx = %E\n', emit);
fprintf(fid, ' emity = %E\n', emit);
fprintf(fid, ' bunch = %E\n', bunch);
fprintf(fid, ' bunchphase = %E\n', bunchphase);
fprintf(fid, ' xlamds = %E\n', lambda);
fprintf(fid, ' prad0 = %E\n', prad0);
fprintf(fid, ' zrayl = %E\n', zrayl);
fprintf(fid, ' zwaist = %E\n', zwaist);
fprintf(fid, ' quadf =  %E\n' , 0);
fprintf(fid, ' quadd =  %E\n', 0);
fprintf(fid, ' fl =  %E\n', 0);
fprintf(fid, ' f1st =  %E\n', 0);
fprintf(fid, ' dl =  %E\n', 0);
fprintf(fid, ' drl =  %E\n', 0);
fprintf(fid, ' nwig = %d\n', nwig);
fprintf(fid, ' ncar = %d\n', ncar);
fprintf(fid, ' dgrid = %E\n', dgrid);
fprintf(fid, ' curpeak = %E\n', curpeak);
fprintf(fid, ' curlen = %E\n', git_par.curlen);
fprintf(fid, ' nbins = 16 \n');
fprintf(fid, ' nscr = %d\n', git_par.nscr);
fprintf(fid, ' nscz = %d\n', git_par.nscz);
fprintf(fid, ' nptr = 20 \n');
fprintf(fid, ' idmpfld = 1 \n');
fprintf(fid, ' idmppar = 1 \n');
fprintf(fid, ' iphsty = 1 \n');
fprintf(fid, ' ishsty = 1 \n');
fprintf(fid, ' ipradi = 0 \n');
fprintf(fid, ' isradi = 0 \n');
fprintf(fid, ' ippart = 0 \n');
fprintf(fid, ' ispart = 0 \n');
fprintf(fid, ' iscan = %d\n', iscan);
fprintf(fid, ' nscan = %d\n', nscan);
fprintf(fid, ' svar = %E\n', svar);
fprintf(fid, ' itdp = %d\n', itdp);
fprintf(fid, ' nslice = %d\n', nslice);
fprintf(fid, ' zsep = %d\n', zsep);
fprintf(fid, ' delz = %d\n', delz);
fprintf(fid, ' iotail = %d\n', iotail);
fprintf(fid, ' itgaus = %d\n', itgaus);
fprintf(fid, ' iorb = 0 \n');
fprintf(fid, ' isravg = 0 \n');
fprintf(fid, ' shotnoise = 1 \n');
fprintf(fid, ' outputfile = "%s"\n', outputfilename);
if (git_par.radfile)
    fprintf(fid, ' radfile = "%s"\n', radfilename);
end
if (quadrupole_lattice)
    fprintf(fid,'%s\n',' magin = 1');
    fprintf(fid,'%s\n',' maginfile = "TESS3.0/mg.mag.in"');
    fprintf(fid,'%s\n',' magoutfile = "TESS3.0/mg.mag.out"');
end
if (git_par.beamfile)
    fprintf(fid, ' beamfile = "%s"\n', beamfilename);
end
if(j>1 & itdp & ~iotail)
    git_par.ntail = git_par.ntail + slippage;
end
fprintf(fid, ' ntail = %d\n', git_par.ntail);
fprintf(fid, '%s\n',' LOUT= 1 0 1 1 1 0 1 1 1 1 1 0 0 1 0 0 0'); 

if (recirculated)
    fprintf(fid, ' fieldfile = "%s.dfl"\n', inputfilename);
end

if j==1 & prebunched
    fprintf(fid, ' fieldfile = "%s.dfl"\n', inputfilename);
    fprintf(fid, ' partfile = "%s.dpa"\n', inputfilename);
    if(itdp)
    offset = git_par.R56slippage;
    end
end
fprintf(fid, ' alignradf = 1 \n');
fprintf(fid, ' offsetradf = %i \n', offset);

if (j > 1)
    fprintf(fid, ' fieldfile = "%s.dfl"\n', inputfilename);
    fprintf(fid, ' partfile = "%s.dpa"\n', inputfilename);
end

fprintf(fid, ' $end\n');
st = fclose('all');

%% Run Genesis using script file in home-directory. Script file avoids recurring call to bash and speeds up the code.

if (git_par.parallel)
system(strcat(dircyg,'bin\bash --login -c "./scriptmpi" > ', dirname,'output.txt'));
else
system(strcat(dircyg,'bin\bash --login -c "./script" > ',dirname,'output.txt'));
end 

output = strcat(dirname,'mg.out');

if ( j > 1 )
clear Cpart Cfield Cfieldt Epart0 C peak ncapt gamma phiabs;
end

j = j + 1;
z(j) = z(j-1) + nwig*lambda_w;
lambda_undulator(j) = lambda_w;
K_undulator(j) = aw0;
fbess = besselj(0,aw0.^2/2/(1+aw0.^2))-besselj(1,aw0.^2/2/(1+aw0.^2));
gamma_res = sqrt(lambda_w/2/lambda*(1+aw0^2));
gammaresonant(j)=gamma_res;

%% slippage and number of slices
if (itdp & ~iotail)
    if (offset)
        slippage = floor(nwig/zsep);
    else 
        slippage = 0;
    end
    slippagetot = slippagetot + slippage;
    niter = niter - slippage;
end

%% Import field data to compute K_laser
outputfieldfile = strcat(dirname,'mg.out.dfl');
Cfieldfile = fopen(outputfieldfile);
Cfieldt = fread(Cfieldfile,2*ncar*ncar*niter,'double');

for js = 1:niter
pmidsc(js) = Cfieldt(2*(ncar*ncar+1)/2+(js-1)*2*ncar*ncar).^2+Cfieldt(2*(ncar*ncar+1)/2+1+(js-1)*2*ncar*ncar).^2;
if (itdp)
p_mid_M(js,j)=pmidsc(js); 
phi_mid_M(js,j)=angle(Cfieldt(2*(ncar*ncar+1)/2+(js-1)*2*ncar*ncar) + 1i* Cfieldt(2*(ncar*ncar+1)/2+1+(js-1)*2*ncar*ncar)); 
end
end
p_mid(j) = max(pmidsc); 
fclose(Cfieldfile);

% Pulse evolution in tapered FEL amplified and oscillator will strongly depend from the choice of target tapering intensity 
% Taper to the max intensity or to the FWHM or to 10%

if(helical) 
    E_laser = sqrt(p_mid(j)/4/dgrid/dgrid*ncar*ncar*377);
else
    E_laser = sqrt(p_mid(j)/4/dgrid/dgrid*ncar*ncar*377*2);
end
K_laser(j) = E_laser/2/pi*lambda/511000;

    
%% Open particle and field file for tapering optimization
outputpartfile = strcat(dirname,'mg.out.dpa');
Cpoint = fopen(outputpartfile);
outputfieldfile = strcat(dirname,'mg.out.dfl');
Cfieldfile = fopen(outputfieldfile);

for islice = 1:niter  % loop on slices
%% read in particleoutput
Cpart = fread(Cpoint,6*npart,'double');

for ip = 1:npart;
    phiabs(ip) = Cpart(ip+npart);
    gamma(ip) = Cpart(ip);
    gammap(ip) = gamma(ip);
    xpart(ip) = Cpart(ip+2*npart);
    ypart(ip) = Cpart(ip+3*npart);
    if correct_R56 
        %phiabs(ip) = phiabs(ip) + 2*pi* (-1+1./gammapart(ip).^2*(1+awd^2)*lambda_w/lambda/2);
        %phiabs(ip) = phiabs(ip) - (1./gammapart(ip).^2)*2*pi*lambda_w/lambda/2;
        phiabs(ip) = phiabs(ip) + (1./gamma(ip).^2- 1./gamma_res.^2)*(awd^2)*2*pi*lambda_w*drift/lambda/2; 
    end
end
energy(slippagetot+islice,j) = mean(gamma);
eta(slippagetot+islice,j)=mean(gamma)/gamma0-1;
bunching(slippagetot+islice,j) = abs(sum(exp(1i.*phiabs))/npart);
bunchphases(slippagetot+islice,j) = angle(sum(exp(i.*phiabs))/npart);
xavg(slippagetot+islice,j) = mean(xpart);
xrms(slippagetot+islice,j) = std(xpart);
yrms(slippagetot+islice,j) = std(ypart);
rpart = sqrt(xpart.^2+ypart.^2);

%% read in fieldoutput
Cfield = fread(Cfieldfile,2*ncar*ncar,'double');
recfield = Cfield(1:2:2*git_par.ncar*git_par.ncar-1);
imcfield = Cfield(2:2:2*git_par.ncar*git_par.ncar);

if git_par.fieldcontourplot ==1 && ((itdp & islice==round(git_par.nslice/2)) || itdp==0)
    field=recfield+1i*imcfield;
    fieldsq=reshape(field,git_par.ncar,git_par.ncar);
    nx = ((1:git_par.ncar)-(git_par.ncar+1)/2);
    dx = 2*git_par.dgrid/git_par.ncar;
    xgridp=nx*dx;
    dx = 2*git_par.dgrid/git_par.ncar;
    fieldcontour=figure(112);
    contourf(xgridp,xgridp,abs(fieldsq));
    cb=colorbar;
    ylabel(cb,'power');
    drawnow
end

I0 = (recfield.^2+imcfield.^2);
power(slippagetot+islice,j) = sum(sum(I0));
Egrid = reshape(I0,ncar,ncar);

phasefield = angle(recfield+1i*imcfield);
phasegrid = reshape(phasefield,ncar,ncar);

peakintensity = max(max(Egrid));
Egrid = Egrid./peakintensity;
gridsize = 2*dgrid/ncar;

% Calculate rms of field profile
%f = fit (xc.',sum(Egrid).','gauss1');
%Egauss = exp(-(xpart./r_gauss(j)).^2/2-(ypart./r_gauss(j)).^2/2);

xc = (1:ncar);
[imax,pmax] = max(sum(Egrid));
r_rms = sqrt(sum((xc-pmax).^2.*sum(Egrid))/sum(sum(Egrid)))*gridsize;
r_size(slippagetot+islice,j) = r_rms;

xgrid = round((xpart+dgrid)/(gridsize));
ygrid = round((ypart+dgrid)/(gridsize));

for ipart = 1:npart
    Epart(ipart) = (Egrid(xgrid(ipart),ygrid(ipart)));
    phasepart(ipart) = phasegrid(xgrid(ipart),ygrid(ipart));
end

if(helical)
    Hfactor = lambda_w/lambda*aw0*K_laser(end);
else
    Hfactor = lambda_w/lambda*aw0*K_laser(end)*fbess/sqrt(2);    
end

theta = mod(phasepart+phiabs-pi,2*pi)-pi;

if(verbose)
    verbosefig=figure(100);
    clf(gcf);
    pointsize=25;
    colorlist=parula(length(gamma));
    scatter(theta,gamma,pointsize,colorlist,'filled')
    xlabel('\phi')
    ylabel('\gamma')
    xlim([-pi pi]);
if git_par.writevid==1
    drawnow;
    getfig=getframe(verbosefig);
    axis([-pi pi 0 gammaresfinal*1.1]) 
    MovieTRIG(j+movieoffset) = getframe(100);
%    writeVideo(vidfile,getfig);
end
end

% First cut
Hamilton = (gamma-gamma_res).^2 - Hfactor.*sqrt(Epart).*(cos(theta)+1);      % calculate the amount of particles trapped if we do not change gamma_res i.e. for psi_res = 0;

csi_res = pi/4;
Hamilton2 = (gamma-gamma_res).^2 - Hfactor.*sqrt(Epart).*(cos(theta)+cos(csi_res)+(theta+csi_res-pi)*sin(csi_res));
sum(Hamilton2<0 & theta<pi-csi_res)/1024

theta(Hamilton>0) = [];
Epart(Hamilton>0) = [];
phasepart(Hamilton>0) = [];
gamma(Hamilton>0) = [];
np = size(gamma,2);

% Second cut. Determine threshold intensity and phase
ncapt(islice) = np;
peak(islice) = peakintensity;
if (np > 0)
    Epart0(islice) = min(Epart)/fac0;      %minimum field felt by particles trapped in the bucket
else
    Epart0(islice) = 0;
end

ncslice(slippagetot+islice,j) = np/npart;
energyslice(slippagetot+islice,j) = mean(gamma);
phislice(slippagetot+islice,j) = phasegrid(ncar/2+0.5,ncar/2+0.5);
thetaslice(slippagetot+islice,j) = mean(theta);
end
fclose(Cpoint);
fclose(Cfieldfile);

%% Define fac and psi_res
psi_res = psi_0;
if(itdp)
%fac = sum(peak.*(sqrt(Epart0)))/sum(peak);                    % should use weighted average. not just mean. 
fac = sum(current.*(sqrt(Epart0)))/sum(current);               % should use weighted average. not just mean. 
np = sum(current.*ncapt)/sum(current);
else
    fac = sqrt(Epart0);
    np = ncapt;
end

resphase(j)=asin(fac*sin(psi_res));
d3fac(j) = fac;

%% Tapering equations. Need fac (intensity fraction) and psi_res (resonant phase)
if(z(end)<z_stop) 


if(period_tapering > 0) % Helical only
    incl = - 8*pi*aw0*K_laser(j)*fac*sin(psi_res)*lambda_w*nwig / (1+aw0^2+lambda_w*2*aw0*dK_und(lambda_w));
    lambda_w = lambda_w+incl;
    aw0 = K_und(lambda_w);
elseif (period_tapering == 0 & ((quadrupole_lattice==1 & lat(1,j-1)) | quadrupole_lattice ==0))
    if(helical)
            incK = -4*pi*sin(psi_res)*K_laser(j)*fac*nwig;
            aw0 = aw0+incK;
    else
            incK = -4*pi*sin(psi_res)*fbess*K_laser(j)*fac*nwig;
            aw0 = aw0+incK/1.414;
    end
elseif (period_tapering == -1)
  aw0 = Nocibur_K_init(j);
  lambda_w = Nocibur_periods(j);
  formatSpec2 = 'Pre-set undulator. aw0 %.3f lambda_w %.3f \n';
end
end


%% Output to command window
trapped(j) = min(ncapt)/npart;
telapsed=toc(tstart);
if mod(j,outputfrequency) == 0 
    formatSpec = '%f out of %f total length. Power %.3e Bunching %.3f resphase %.3f fac %.3f trap %.3f in %.3f sec \n';
    fprintf(formatSpec,z(j),z_stop,max(power(max(1,slippagetot):end,end)),max(bunching(max(1,slippagetot):end,end)),psi_res,fac,np/npart,telapsed*outputfrequency);
end

%% Correct theta and/or correct R_56 in particle files
if ( j >= 2)
    Cpoint = fopen(outputpartfile,'r+');
    Ctot = fread(Cpoint,6*npart*niter,'double');
    fclose(Cpoint);
    avgpsi_niter=zeros(1,niter);
        for islice = 1:niter
            for ip = 1:npart
                    phiabs(ip) = Ctot(ip+npart+(islice-1)*6*npart);
                    gammapart(ip) = Ctot(ip+(islice-1)*6*npart); 
                    if correct_R56
                        %phiabs(ip) = phiabs(ip) + 2*pi* (-1+1./gammapart(ip).^2*(1+awd^2)*lambda_w*drift/lambda/2);
                        %phiabs(ip) = phiabs(ip) - (1./gammapart(ip).^2)*2*pi*lambda_w*drift/lambda/2;
                        phiabs(ip) = phiabs(ip) + (1./gammapart(ip).^2- 1./gamma_res.^2)*(awd^2)*2*pi*lambda_w*drift/lambda/2; 
                    end
                    if correct_theta
                    phiabs(ip) = phiabs(ip)+deltatheta(j);
                    end
                    Ctot(ip+npart+(islice-1)*6*npart) = phiabs(ip);                    
            end
            avgpsi_niter(islice)=mean(phiabs);
        end
    avgpsi(j)=mean(avgpsi_niter);      
    Cpoint = fopen(outputpartfile,'w+');
    Ccount = fwrite(Cpoint,Ctot,'double');
    fclose(Cpoint);
end

if git_par.plot_phase==1 & git_par.quadrupole_lattice==1&&(at_undexit==1 && prebunched==1)  %plot at the entrance of prebuncher and after each undulator exit
plot_part_phase;
end

end


 


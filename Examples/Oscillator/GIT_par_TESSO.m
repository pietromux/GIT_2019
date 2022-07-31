%% README

% In order to minimize the number of system calls we use script files to
% run genesis from Matlab

% in homedir the following scripts file must be created
% scriptmpi
% >mv TESS3.0/mg.out.dfl TESS3.0/mg.in.dfl
% >mv TESS3.0/mg.out.dpa TESS3.0/mg.in.dpa
% >sleep 0.05 
% >mpirun -verbose -np 7 TESS3.0/genesis_mpi.exe < TESS3.0/INPUT.txt
% 
% script
% >mv TESS3.0/mg.out.dfl TESS3.0/mg.in.dfl
% >mv TESS3.0/mg.out.dpa TESS3.0/mg.in.dpa
% >sleep 0.05 
% >genesis.exe < TESS3.0/INPUT.txt

% in dirname the INPUT.txt file must be created
% >TESS3.0/mg.in

clear all 
close all

% GIT_dir should be in MATLAB path (or in same directory as input file)
GITMATLAB_path= genpath('C:\Users\pietro\Desktop\GIT_2019');
addpath(GITMATLAB_path);
GIT_dir;

git_par.itdp = 0;                %1: time-dependent , 0: time-independent
    
git_par.f = 3.47;               % focal length of cavity lens (m)
git_par.d = .33;                % distance of the lens from the center of the cavity (m)

saveplots=0;               %saves TESSO plots in figdir
% if git_par.itdp==1
% figdirname=sprintf('cd=%.f trans=%.3f Q=%.2e f=%.3f d=%.3f',git_par.cavitydetuning,git_par.transmission,git_par.Q,git_par.f,git_par.d);
% else
% figdirname=sprintf('trans=%.3f',git_par.transmission);
% end

if saveplots==1
figdir=[datadir,'TESSOfigures\',datetimenow,'_',figdirname,'\'];
mkdir(figdir);
end
%% Parameter list to run GIT_fun 
npasses = 20;
for iss = 1:npasses
    tstart=tic;
    if iss==1
    recirculated = 0;
    git_par.radfile = 1;
    else
    recirculated = 1; 
    git_par.radfile = 0;
    fprintf('recirc %d pass # %d \n', recirculated, iss);
    end
    
%% Parameter list to run GIT_fun 
git_par.prad0 = 1e9;    
git_par.parallel = 0;            % 1 parallel ; 0 single processor 
git_par.helical = 1;             % 0 for planar; 1 for helical
git_par.period_tapering = 0;    % 0 for constant period, 1 for varying period, -1 for predetermined period
git_par.verbose = 0;             % plot phase space during run
git_par.outputfrequency = 1;     % # of steps between command window status report 
git_par.plot_phase=0;
git_par.fieldcontourplot=1;        %output central contour plot after each Genesis run 
%% Different TESSA optimization schemes
git_par.lps_predictor = 0;       % predicts phase space to decide tapering step. only works for NWIG small - obsolete
git_par.correct_theta = 0;       % every drift section a phase shifter pushes the trapped phasespace in the center of the resonant bucket
git_par.correct_R56 =   1;       % Correct for artificial r56 introduced by Genesis in drifts 

git_par.psi_0 = 0.9;               % resonant phase
git_par.fac0 = 1;                   % For uniform distribution use fac0 = 1 psi_0 = 0.9 
                                            % For parabolic distribution use fac0 = 0.95 psi_0 = 1.3      

git_par.offset = 0;                      % slippage control parameter. should be -1 - unless this is a full undulator simulation
 
git_par.quadrupole_lattice = 2;   % 1 Use quadrupole lattice from magin file
                                                   % 0 Keep a constant beam size using xkx and xky

git_par.drift = 25;                  % drift length in undulator period units 
git_par.rxbeam = 30.00e-6;   % Matched beam radius
git_par.rybeam = 30.00e-6;   % Matched beam radius

%% GENESIS parameters
git_par.itgaus = 2;     % 0 parabolic ; 1 gaussian ; 2 uniform transverse profiles
git_par.iscan = 0;
git_par.nscan = 1;
git_par.svar = 0.00001;
git_par.delz = 1;
git_par.zsep = 1;       % zsep = nwig for itdp = 1
git_par.iotail = 1;
git_par.slippage = 0;
git_par.slippagetot = 0;

git_par.wcoefz = [ 0 0 0 ];
git_par.lambda_w0 = 0.032;
git_par.npart = 256;
git_par.gamma0 = 733.8;  
git_par.delgamma = git_par.gamma0*0.001;
git_par.emit = 2e-6;
git_par.ntail = 0;
git_par.bunch = 0.0;
git_par.bunchphase = -0; 
git_par.lambda = 0.266e-6;
git_par.zrayl =  1.3508;
git_par.zwaist = 1.4844;
git_par.nwig = 1;
git_par.ncar = 131;
git_par.dgrid = 3.2e-3;
git_par.nscz = 0;
git_par.nscr = 0;
git_par.pulselength = 180e-6;    % Laser pulse length
git_par.beamfile = 1;

%% parameters used for period tapering 
git_par.gap = 0.003;
git_par.Brem = 1.5;

%% parameters for transverse matching. Calculated separately.   
%%%48T/m 227T/m parameter

git_par.rxbeam_match=94.95e-6;
git_par.rybeam_match=115.23e-6;

git_par.alphax= -0.822;
git_par.alphay=  1.577;
git_par.quadgradient = 0;
git_par.xkx = 0.5;
git_par.xky = 0.5;
git_par.lquad = 0;
git_par.fododrift = 0;

%% Run GIT_fun 

%% Tapered undulator
clear z j;
if (git_par.period_tapering < 0)
temp = load('outres.mat');
Nocibur_periods = circshift(temp.lambda_undulator,-1,2);
Nocibur_K_init = circshift(temp.K_undulator,-1,2);
Nocibur_periods(end) = [];
Nocibur_K_init(end) = [];
end

git_par.nslice=800;
git_par.cavitydetuning = 170;    % In units of zsep
git_par.transmission = .2;          % Power transmission through one cavity pass
git_par.sigmaomega = 0.005*git_par.nslice*git_par.zsep;
git_par.Q= 1/ git_par.sigmaomega;
git_par.curpeak = 1e3;
git_par.curlen = 90e-6;     % Electron beam pulse length
curpeakpos=git_par.nslice/2;
radshift=curpeakpos-120;
git_par.charge=300e-12;
write_radfile;
p_rad_0=p_rad;
beamtype=0;
write_beamfile;

%% prebuncher
git_par.nwig=53;
git_par.quadrupole_lattice=2;
git_par.z_stop =.032*53;          % Simulation length
prebunched = 0;
git_par.period_tapering = 0;
git_par.psi_0 = 0;            
totalmagfilename=['mg_pb.mag.in'];
copyfile(totalmagfilename,[dirname,'mg.mag.in']);
git_par.awd=0;
git_par.pmidvar=1;
readin_totalmagfile;
git_par.drift = 45;                  % drift length in undulator period units 
GIT_fun_2019;
p_rad_1=power(:,end);
%output_onerun_plots=1;
%GIT_fun_output_onerun;
%saveas(gcf,[figdir,'poutPB_',num2str(iss),'.png']);
% figure(pmidplot)
% saveas(gcf,[figdir,'poutPBpmid_',num2str(iss),'.png']);
%figure(fieldcontour)
%saveas(gcf,[figdir,'poutPBfieldcontour_',num2str(iss),'.png']);

init_rms(iss)=r_rms; %plot initial rad spot size

R56buncher = 10e-6;
phaseshift = 1;
prop_part;
plot_part_phase;
% undulator
git_par.correct_theta = 1;       % every drift section a phase shifter pushes the trapped phasespace in the center of the resonant bucket
git_par.correct_R56 =   1;       % Correct for artificial r56 introduced by Genesis in drifts 

% clearvars -except git_par phasej num radshift niter R56buncher totaltimestart pbdir p_rad p_rad_0 p_rad_1 energy;
clear power energy bunching j z 
prebunched = 1;
git_par.period_tapering = 0;
git_par.z_stop = 1.152;         
git_par.R56slippage = -2*round(R56buncher/git_par.lambda);
radshift=radshift-git_par.R56slippage;
write_radfile;
git_par.prad0 = 1e9;
git_par.nwig=36;

%name of the four magin files for each undulator section containing nominal
%tapering
maginfile1=['mg_und_1.mag.in'];
maginfile2=['mg_und_2.mag.in'];
maginfile3=['mg_und_3.mag.in'];
maginfile4=['mg_und_4.mag.in'];

git_par.drift = 7;                  % drift length in undulator period units 

dtheta1=-.7166;
dtheta2=-.6916;
dtheta3=-.2479;
figure(113)
subplot(3,2,2)
contourf(xgridp,xgridp,abs(fieldsq));
legend('buncher');


%%first undulator
git_par.deltatheta=[0,dtheta1];
copyfile(maginfile1,[dirname,'mg.mag.in']);
GIT_fun_2019;
git_par.R56slippage=-1; %slippage after the first undulator
%output_onerun_plots=0; %do not output plots for the single runs and just read the parameter
%read_genesis_onerun;
figure(113)
subplot(3,2,3)
contourf(xgridp,xgridp,abs(fieldsq));
legend('1st und');

power1=power;z1=z;energy1=energy;bunching1=bunching; r_size1=r_size; K_undulator1=K_undulator; xrms1=xrms; yrms1=yrms; 
%%second undulator
copyfile(maginfile2,[dirname,'mg.mag.in']);
git_par.deltatheta=[0,dtheta2];
GIT_fun_2019;
%read_genesis_onerun;
figure(113)
subplot(3,2,4)
contourf(xgridp,xgridp,abs(fieldsq));
legend('2nd und');

power2=power;z2=z;energy2=energy;bunching2=bunching; r_size2=r_size; K_undulator2=K_undulator; xrms2=xrms; yrms2=yrms; 
%%third undulator
copyfile(maginfile3,[dirname,'mg.mag.in']);
git_par.deltatheta=[0,dtheta3];
GIT_fun_2019;
%read_genesis_onerun;
figure(113)
subplot(3,2,5)
contourf(xgridp,xgridp,abs(fieldsq));
legend('3rd und');

power3=power;z3=z;energy3=energy;bunching3=bunching; r_size3=r_size; K_undulator3=K_undulator; xrms3=xrms; yrms3=yrms; 
%%fourth undulator
copyfile(maginfile4,[dirname,'mg.mag.in']);
GIT_fun_2019;
%read_genesis_onerun;
figure(113)
subplot(3,2,6)
contourf(xgridp,xgridp,abs(fieldsq));
legend('4th und');

power4=power;z4=z;energy4=energy;bunching4=bunching; r_size4=r_size; K_undulator4=K_undulator; xrms4=xrms; yrms4=yrms; 
%%put together outputs from the four simulation
power=[power1(:,2:end),power2(:,2:end),power3(:,2:end),power4(:,2:end)];
z=[z1(2:end);z2(2:end)+1.152;z3(2:end)+1.152*2;z4(2:end)+1.152*3];
energy=[energy1(:,2:end),energy2(:,2:end),energy3(:,2:end),energy4(:,2:end)];
bunching=[bunching1(:,2:end),bunching2(:,2:end),bunching3(:,2:end),bunching4(:,2:end)];
r_size=[r_size1(:,2:end),r_size2(:,2:end),r_size3(:,2:end),r_size4(:,2:end)];
K_undulator=[K_undulator1(2:end);K_undulator2(2:end);K_undulator3(2:end);K_undulator4(2:end)];
xrms=[xrms1(:,2:end),xrms2(:,2:end),xrms3(:,2:end),xrms4(:,2:end)];
yrms=[yrms1(:,2:end),yrms2(:,2:end),yrms3(:,2:end),yrms4(:,2:end)];

% if itdp==0
% GIT_fun_output_tindep;
% subplot(2,3,4)
% ylim([min(K_undulator(K_undulator~=0)) max(K_undulator)])
% else
% GIT_fun_output_tdep_sub;
% subplot(3,3,4)
% xlim([266e-9-5e-9 266e-9+5e-9])
% end
      
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
% outres.xrms = xrms;
% outres.d3fac = d3fac;
% outres.resphase = resphase;
% outres.gammaresonant = gammaresonant;
% outres.theta = thetaslice;

%energy efficiency -- time-indep Erad calculation assumes rectangular
%current of given charge and curpeak

%% TESSO output ;

Ebeam=git_par.gamma0 *git_par.charge*.511e6;
Eloss=(energy(end)-energy(1))*git_par.charge*.511e6;
Erad=(sum(power(:,end))-sum(power(:,1)))*git_par.lambda*git_par.zsep/3e8;
Ueff = Erad/Ebeam

rad_vs_und(:,iss)=mean(power,1);
rad_vs_beam(:,iss)=power(:,end);
blist(:,iss)=bunching(:,end);
Eff(iss)=Ueff;
if(~itdp)
    Eff(iss) = (sum(power(:,end))-sum(power(:,1)))/(git_par.gamma0*.511e6*git_par.curpeak);
end

%transmission and filtering
write_recirc_field;

%optical setup for focusing the beam
if (iss ~= npasses)
    prop_spectr_recirc;
end
figure(113)
subplot(3,2,1)
contourf(abs(reshape(propfield_m,ncar,ncar)))
legend('input')

if saveplots==1
saveas(outfig,[figdir,'poutUndExit_',num2str(iss),'.png']);
saveas(fieldcontour,[figdir,'poutUndExitfieldcontour_',num2str(iss),'.png']);
end
%copyfile([dirname,'mg.out.dfl'],[dirname,'TESSO_tdep_rec_',num2str(iss),'.out.dfl']); %after above script uses mg.out.dfl to run the prebuncher replace back with the saved dfl
toc(tstart)
end
rmsfig=figure(14);
plot(init_rms);
title('initial radiation spot size at prebuncher entrance')
xlabel('iteration')
ylabel('r-size(m)')

maxpowerfig=figure(15);
plot(rad_vs_beam);
xlabel('coordinate along the bunch')
ylabel('max power(W)')
title('power at the undulator exit')

bunchfig=figure(16);
plot(blist);
title('bunching factor at the undulator exit')
xlabel('coordinate along the bunch')
ylabel('bunching factor')

efffit = figure(17)
plot(Eff);
title('efficiency vs. number of passes')
xlabel('iteration')
ylabel('efficiency')


if saveplots==1
saveas(rmsfig,[figdir,'initrms.png']);
saveas(maxpowerfig,[figdir,'pmax.png']);
saveas(bunchfig,[figdir,'blist.png']);
end

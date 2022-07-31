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

%% Parameter list to run GIT_fun 
git_par.prad0 = 1e2;    
git_par.parallel = 0;            % 1 parallel ; 0 single processor 
git_par.helical = 0;             % 0 for planar; 1 for helical
git_par.period_tapering = -1;    % 0 for constant period, 1 for varying period, -1 for predetermined period

git_par.verbose = 0;             % plot phase space during run
git_par.outputfrequency = 1;     % # of steps between command window status report 
git_par.fieldcontourplot = 0;    % plot central slice transverse profile contour plot at the end of each Genesis run
git_par.writevid=0;                     % write a video file with phasespace evolution

%% Different TESSA optimization schemes
git_par.lps_predictor = 0;       % predicts phase space to decide tapering step. only works for NWIG small - obsolete
git_par.correct_theta = 0;       % every drift section a phase shifter pushes the trapped phasespace in the center of the resonant bucket
git_par.correct_R56 =   0;       % Correct for artificial r56 introduced by Genesis in drifts 

git_par.psi_0 = 0.0;               % resonant phase
git_par.fac0 = 1.0;                   % For uniform distribution use fac0 = 1 psi_0 = 0.9 
                                            % For parabolic distribution use fac0 = 0.95 psi_0 = 1.3      

git_par.offset = 0;                      % slippage control parameter. should be -1.
git_par.quadrupole_lattice = 0;   % 1 Use quadrupole lattice from magin file. Simulate period by period to optimize
                                                   % 0 Keep a constant beam size using xkx and xky
                                                   % 2 Simulate long undulator sections ending with drift (no optimization)

git_par.drift = 25;                    % drift length in undulator period units 
git_par.rxbeam = 10.00e-6;   % Matched beam radius (quadrupole_lattice = 0)
git_par.rybeam = 10.00e-6;   % Matched beam radius

%% GENESIS parameters
git_par.itgaus = 1;     % 0 parabolic ; 1 gaussian ; 2 uniform transverse profiles
git_par.iscan = 0;
git_par.nscan = 1;
git_par.svar = 0.1;
git_par.delz = 1;
git_par.itdp = 0;
git_par.zsep = 1;       % zsep = nwig for itdp = 1
git_par.iotail = 1;
git_par.slippage = 0;
git_par.slippagetot = 0;

git_par.wcoefz = [ 0 0 0 ];
git_par.lambda_w0 = 0.032;
git_par.npart = 1024;
git_par.gamma0 = 15.656;  
git_par.delgamma = git_par.gamma0*1e-6;
git_par.emit = 0.01e-6;
git_par.curpeak = 60;
git_par.curlen = 90e-6;     % Electron beam pulse length
git_par.ntail = 0;
git_par.bunch = 0.0;
git_par.bunchphase = 0; 
git_par.lambda = 7.02e-4;
git_par.zrayl =  0.0108;
git_par.zwaist = 0.00;
git_par.nwig = 1;
git_par.ncar = 211;
git_par.dgrid = 1.185e-3*sqrt(2);
%git_par.dgrid = 1.596e-3;
git_par.nscz = 0;
git_par.nscr = 0;
git_par.pulselength = 180e-6;    % Laser pulse length
git_par.beamfile = 1;
git_par.radfile = 1;

%% Parameters used for period tapering 
%%%%% Ignored if period_tapering < 1 %%%%%
git_par.gap = 0.003;
git_par.Brem = 1.5;

%% parameters for transverse matching. Calculated separately. 
%%%%%%%% Ignored if quadrupole_lattice == 0 %%%%%%%%
git_par.rxbeam_match=94.9e-06;
git_par.rybeam_match=115.23e-06;

git_par.alphax=-0.822;
git_par.alphay=1.577;

git_par.quadgradient = 0;
git_par.xkx = 0.5;
git_par.xky = 0.5;
git_par.lquad = 0;
git_par.fododrift = 0;
git_par.plot_phase=0;

%% Time-dependent gaussian/rectangular beam
git_par.nslice=400;
curpeakpos=git_par.nslice/2;
radshift=curpeakpos-120;
write_radfile;
git_par.charge=300e-12;
p_rad_0=p_rad;
beamtype=0;
write_beamfile;
if git_par.writevid==1
 vidfile = VideoWriter([datadir,'phasemovie.mp4'],'MPEG-4');
 open(vidfile);
end
%% THz FEL run GIT_fun
recirculated = 0;
prebunched = 0;
git_par.dgrid = 1.6e-3;
git_par.lambda = 6e-4;
nd = 11;
nlam = 21;
nperiods = 60;
dgrid0 = git_par.dgrid;
lambda0 = git_par.lambda;
for j1 = 1:nd
 for kk = 1:nlam
    tic
   clear power bunching
    git_par.dgrid = dgrid0+(j1-(nd+1)/2)*5e-5;   
    git_par.lambda = lambda0 + 0.2e-4*(kk-(nlam+1)/2);
    fprintf('%f %f',git_par.dgrid,git_par.lambda);
    git_par.aw0 = 2.09155;
    git_par.awd = git_par.aw0;
    git_par.z_stop =git_par.lambda_w0*nperiods;          % Simulation length
    git_par.nwig = nperiods;
    git_par.drift = 0; 
    np = git_par.z_stop / (git_par.lambda_w0*git_par.nwig);
    Nocibur_periods(1:np+1) = git_par.lambda_w0;
    Nocibur_K_init(1:np+1) = git_par.aw0;
    GIT_fun_2019;
    nheader = 187+nperiods;
    output = 'C:\cygwin\home\pietro\TESS3.0\mg.out';
    clear power z r_size
    C = importdata(output,' ',nheader);
    power = C.data(:,1);
    r_size = C.data(:,5);
    z = (0:nperiods)*git_par.lambda_w0;
    res(j1,kk) = power(end);
    [res3(j1,kk),res4(j1,kk)] = max(power);
    res5(j1,kk) = r_size(res4(j1,kk));
    res6(j1,kk) = std(gammapart);
    la(kk) =  git_par.lambda;
    gr(j1) = git_par.dgrid;
    toc(tic)
    fitregion = z>0.3 & z<z(res4(j1,kk));
    B = polyfit(z(fitregion),log(power(fitregion)), 1);
    PowerLgain = 1/B(1);
    res2(j1,kk) = PowerLgain;
    figure(1)
    semilogy(z,power);
    hold on
    semilogy(z,exp(polyval(B,z)))
 end 
end
%%  Resonance curve %%
figure(2)
[mgla,mggr] = meshgrid(la,gr);
surf(mgla,mggr,res)
%GIT_fun_output;
figure(3)
surf(mgla,mggr,1./res2)
%GIT_fun_output;
figure(4)
surf(mgla,mggr,res3)
figure(5)
surf(mgla,mggr,res4)
figure(6)
surf(mgla,mggr,res5)
%%
out.mgla = mgla;
out.mggr = mggr;
out.res = res;
out.res2 = res2;
out.res3 = res3;
out.res4 = res4;
out.res5 = res5;
save('out')



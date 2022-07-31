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
git_par.prad0 = 1e9;    
git_par.parallel = 0;            % 1 parallel ; 0 single processor 
git_par.helical = 1;             % 0 for planar; 1 for helical

git_par.verbose = 0;                   % plot phase space during run
git_par.outputfrequency = 1;     % # of steps between command window status report 
git_par.fieldcontourplot = 0;      % plot central slice transverse profile contour plot at the end of each Genesis run
git_par.writevid=0;                     % write a video file with phasespace evolution

%% Different TESSA optimization schemes
git_par.lps_predictor = 0;       % predicts phase space to decide tapering step. only works for NWIG small - obsolete
git_par.correct_theta = 0;       % every drift section a phase shifter pushes the trapped phasespace in the center of the resonant bucket 
git_par.correct_R56 =   1;       % Correct for artificial r56 introduced by Genesis in drifts 

git_par.psi_0 = 1.3;               % resonant phase
git_par.fac0 = 1.0;                   % For uniform distribution use fac0 = 1 psi_0 = 0.9 
                                            % For parabolic distribution use fac0 = 0.95 psi_0 = 1.3      

git_par.offset = -1;                      % slippage control parameter. should be -1.
git_par.quadrupole_lattice = 1;   % 1 Use quadrupole lattice from magin file. Simulate period by period to optimize
                                                   % 0 Keep a constant beam size using xkx and xky
                                                   % 2 Simulate long undulator sections ending with drift (no optimization)

%% GENESIS parameters
git_par.itgaus = 1;     % 0 parabolic ; 1 gaussian ; 2 uniform transverse profiles
git_par.iscan = 0;
git_par.nscan = 1;
git_par.svar = 0.00001;
git_par.delz = 1;
git_par.itdp = 0;
git_par.zsep = 1;       % zsep = nwig for itdp = 1
git_par.iotail = 1;
git_par.slippage = 0;
git_par.slippagetot = 0;

git_par.wcoefz = [ 0 0 0 ];
git_par.lambda_w0 = 0.032;
git_par.npart = 1024;
git_par.gamma0 = 733.8;  
git_par.delgamma = git_par.gamma0*0.001;
git_par.emit = 2e-6;
git_par.curpeak = 1000;
git_par.curlen = 90e-6;     % Electron beam pulse length
git_par.ntail = 0;
git_par.bunch = 0.6;
git_par.bunchphase = 0.9; 
git_par.lambda = 0.266e-6;
git_par.zrayl =  1.4;
git_par.zwaist = 1.4;
git_par.nwig = 1;
git_par.ncar = 151;
git_par.dgrid = 3e-3;
git_par.nscz = 0;
git_par.nscr = 0;
git_par.pulselength = 180e-6;    % Laser pulse length
git_par.beamfile = 1;
git_par.radfile = 1;

%% parameters used for period tapering 
git_par.gap = 0.00721;
git_par.Brem = 1.2;

%% parameters for transverse matching. Calculated separately.   
%%%48T/m 227T/m parameter
%c;
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

%% Run GIT_fun 
recirculated = 0;

%% Gaussian/rectangular beam
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

%% Tapered undulator
prebunched = 0;
git_par.drift = 0;
git_par.awd = 0;
git_par.nwig=1;
git_par.quadrupole_lattice = 0;  
git_par.rxbeam = 40.40e-6;   % Matched beam radius (quadrupole_lattice = 0)
git_par.rybeam = 40.40e-6;   % Matched beam radius

git_par.period_tapering = 0;
git_par.z_stop =.032*29*4;         
GIT_fun_2019;
GIT_fun_output;

if itdp==1
    GIT_fun_output_tdep_sub;
else
    GIT_fun_output_tindep;
end

if git_par.writevid==1
    close(vidfile);
end

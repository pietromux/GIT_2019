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

%% GIT_dir should be in MATLAB path (or in same directory as input file)
GITMATLAB_path= genpath('C:\Users\pietro\Desktop\GIT_2019');
addpath(GITMATLAB_path);
GIT_dir;

%% Parameter list to run GIT_fun 
git_par.prad0 = 1;    
git_par.parallel = 0;                 % 1 parallel ; 0 single processor 
git_par.helical = 0;                   % 0 for planar; 1 for helical
git_par.period_tapering = 0;    % 0 for constant period, 1 for varying period, -1 for predetermined period
git_par.verbose = 1;                 % plot phase space during run
git_par.outputfrequency = 1;   % # of steps between command window status report 

% Different TESSA optimization schemes
git_par.lps_predictor = 0;       % predicts phase space to decide tapering step. only works for NWIG small
git_par.correct_theta = 0;       % every drift section a phase shifter pushes the trapped phasespace in the center of the resonant bucket
git_par.correct_R56 =   0;       % Correct for artificial r56 introduced by Genesis in drifts 

git_par.psi_0 = 0;                 % resonant phase
git_par.fac0 = 1;                   % For uniform distribution use fac0 = 1 psi_0 = 0.9 
                                            % For parabolic distribution use fac0 = 0.95 psi_0 = 1.3      

git_par.rxbeam = 66.00e-6;   % Matched beam radius
git_par.rybeam = 66.00e-6;   % Matched beam radius

git_par.offset = 0;                         % slippage control parameter. should be -1.
 
git_par.quadrupole_lattice = 0;   % 1 Use quadrupole lattice from magin file
                                                   % 0 Keep a constant beam size using xkx and xky
                                                   % 2 Use long undulator sections
                                                   
git_par.dircyg = 'C:\cygwin\';
git_par.homedir = strcat(git_par.dircyg,'home\pietro\');
git_par.dirname = strcat(git_par.homedir,'TESS3.0\');

git_par.z_stop = 4.0;          % Simulation length

%% GENESIS parameters
git_par.itgaus = 1;     % 0 parabolic ; 1 gaussian ; 2 uniform transverse profiles
git_par.iscan = 0;
git_par.nscan = 1;
git_par.svar = 0.00001;
git_par.delz = 1;
git_par.itdp = 0;
git_par.iotail = 0;
git_par.nslice = 80;
git_par.zsep = 1;       % It has to be commensurate with nwig: zsep = integer * nwig
git_par.slippage = 0;
git_par.slippagetot = 0;
git_par.wcoefz = [ 0 0 0 ];
git_par.lambda_w0 = 0.018;
git_par.npart = 4096;
git_par.gamma0 = 140;  
git_par.delgamma = git_par.gamma0*0.0001;
git_par.emit = 2e-6;
git_par.curpeak = 200;
git_par.curlen = 0.3e-3;     % Electron beam pulse length
git_par.ntail = 0;
git_par.bunch = 0.0;
git_par.bunchphase = 0; 
git_par.lambda = 0.857e-6;
git_par.zrayl = 0.2;
git_par.zwaist = 0.2;
git_par.nwig = 10;              % This should be 1 unless quadrupole_lattice == 0
git_par.ncar = 191;
git_par.dgrid = 3e-3;
git_par.nscz = 0;
git_par.nscr = 0;
git_par.pulselength = 1e-3;    % Laser pulse length
git_par.beamfile = 0;
git_par.radfile = 0;


%% Run GIT_fun (should be in same directory or in Matlab path)
recirculated = 0;
git_par.plot_phase=1;
git_par.fieldcontourplot = 1;      % plot central slice transverse profile contour plot at the end of each Genesis run
git_par.writevid=0;                     % write a video file with phasespace evolution
prebunched = 0;
zwaist = git_par.zwaist;
git_par.awd=0;

git_par.nwig = 1;
git_par.z_stop = 0.015;         % Simulation length
R56buncher = 0;
phaseshift  = 0;
GIT_fun_2019;

%prop_part_conditioning;
git_par.nwig = 2;
git_par.z_stop = 6;               % Simulation length
prebunched = 1;
git_par.R56slippage = 0;
GIT_fun_2019;
GIT_fun_output;

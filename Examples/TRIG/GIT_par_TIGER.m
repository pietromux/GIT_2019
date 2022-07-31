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

git_par.prad0 = 10e12;    
git_par.parallel = 0;            % 1 parallel ; 0 single processor 
git_par.helical = 1;             % 0 for planar; 1 for helical
git_par.period_tapering = -1;    % 0 for constant period, 1 for varying period, -1 for predetermined period
git_par.verbose = 1;             % plot phase space during run
git_par.outputfrequency = 1;     % # of steps between command window status report 
git_par.fieldcontourplot = 0;      % plot central slice transverse profile contour plot at the end of each Genesis run
git_par.writevid=0;                     % write a video file with phasespace evolution
git_par.plot_phase = 0;     

% Different TESSA optimization schemes
git_par.lps_predictor = 0;       % predicts phase space to decide tapering step. only works for NWIG small
git_par.phase_locking = 0;     % Not benchmarked yet
git_par.correct_theta = 0;       % every drift section a phase shifter pushes the trapped phasespace in the center of the resonant bucket
git_par.correct_R56 =   0;       % Correct for artificial r56 introduced by Genesis in drifts 
git_par.AD_correction = 1;
git_par.change_psires = 0;
git_par.deltapsi = 0.0;         % Dephasing for drifts. Possibly redundant if correct_theta = 1;

git_par.psi_0 = -0.85;            % resonant phase
git_par.fac0 = 1;                   % For uniform distribution use fac0 = 1 psi_0 = 0.9 
                                            % For parabolic distribution use fac0 = 0.95 psi_0 = 1.3      

                                            git_par.rxbeam = 100.00e-6;   % Matched beam radius
git_par.rybeam = 100.00e-6;   % Matched beam radius

git_par.offset = 0;                         % slippage control parameter. should be -1.
 
git_par.quadrupole_lattice = 0;   % Use quadrupole lattice

%% GENESIS parameters

git_par.itgaus = 1;     % 0 pfarabolic ; 1 gaussian ; 2 uniform transverse profiles
git_par.iscan = 0;
git_par.nscan = 1;
git_par.svar = 0.00001;
git_par.delz = 0.01;
git_par.itdp = 0;
git_par.iotail = 0;
git_par.nslice = 80;
git_par.zsep = 1;       % It has to be commensurate with nwig: zsep = integer * nwig
git_par.slippage = 0;
git_par.slippagetot = 0;
git_par.wcoefz = [ 0 0 0 ];
git_par.lambda_w0 = 0.0215;
git_par.npart = 256;
git_par.gamma0 = 197;  
git_par.delgamma = git_par.gamma0*0.0004;
git_par.emit = 1e-7;
git_par.curpeak = 1000;
git_par.curlen = 0.3e-3;     % Electron beam pulse length
git_par.ntail = 0;
git_par.bunch = 0.8;
git_par.bunchphase = -1.05; 
%git_par.prad0 = 30e9;
git_par.lambda = 400e-9;
git_par.zrayl = 0.625;
git_par.zwaist = 1.02;
git_par.nwig = 1;
git_par.ncar = 151;
git_par.dgrid = 4e-3;
git_par.nscz = 0;
git_par.nscr = 0;
git_par.pulselength = 1e-3;    % Laser pulse length
git_par.beamfile = 0;
git_par.radfile = 0;

%% parameters used for period tapering 
git_par.gap = 0.003;
git_par.Brem = 1.5;

%% parameters for transverse matching. Calculated separately.   
git_par.rxbeam_match = 300e-6;
git_par.rybeam_match = 300e-6;
git_par.alphax = 3;
git_par.alphay = 3;
git_par.quadgradient = 0;
git_par.xkx = 0.5;
git_par.xky = 0.5;
git_par.lquad = 0;
git_par.fododrift = 0;

%% Run GIT_fun (should be in same directory or in Matlab path)
recirculated = 0;
%% Prebuncher section
% prebunched =0;
% git_par.period_tapering = -1;
% git_par.nwig = 2;
% R56buncher = 0.00006;   % drift R56 + chicane R56
% phaseshift = -1.2;
% drift_buncher = 0.26;
% bunchernoperiods = 1;
% %write_beamfile;
% %write_radfile;
% git_par.zwaist = git_par.zwaist + 1.0;
% git_par.z_stop = 0.12;
% Nocibur_periods = [0.12, 0.12];
% Nocibur_K_init = [2.52/sqrt(2), 2.52/sqrt(2)];
% M = [1 drift_buncher ; 0 1];
% Mag = 1;
% factor_gamma = 0.25;
% GIT_fun_2016;
% prop_spectr;
% prop_part;
% %% Second buncher
% prebunched = 1;
% factor_gamma = 1;
% git_par.nwig = 3;
% git_par.nslice = git_par.nslice-1;
% %clear z j;
% R56buncher = 0.00002;    % drift R56 + chicane R56
% phaseshift = 1.5;
% drift_buncher = 0.40;
% git_par.z_stop = 0.05;
% Nocibur_periods = [0.05, 0.05];
% Nocibur_K_init = [4.243/sqrt(2), 4.243/sqrt(2)];
% M = [1 drift_buncher; 0 1];
% Mag = 1;
% GIT_fun_2016;
% prop_spectr;
% prop_part;

%% Tapered undulator
prebunched = 0;
git_par.awd=0;
git_par.period_tapering = -1;
zwaist = git_par.zwaist;
delz = 0.01;
clear z j;
if (git_par.period_tapering < 0)
temp = load('-ascii','outres.txt');
Nocibur_periods = circshift(temp(:,1),-1,2)';
Nocibur_K_init = circshift(temp(:,2),-1,2)';
Nocibur_periods(end) = [];
Nocibur_K_init(end) = [];
Nocibur_periods(end) = [];
Nocibur_K_init(end) = [];
secondsect_periods = flip(Nocibur_periods);
secondsect_K = flip(Nocibur_K_init);
gammaresfinal = sqrt(Nocibur_periods(end)*(1+Nocibur_K_init(end)^2)/2/git_par.lambda);
Kdephase = sqrt(1.6*git_par.lambda*gammaresfinal^2/Nocibur_periods(end));
Nocibur_periods = [Nocibur_periods Nocibur_periods(end) secondsect_periods];
Nocibur_K_init = [Nocibur_K_init Kdephase secondsect_K];
Nocibur_periods(end+1) = Nocibur_periods(end);
Nocibur_K_init(end+1) = Nocibur_K_init(end)*0.9;
Nocibur_periods(end+1) = Nocibur_periods(end);
Nocibur_K_init(end+1) = Nocibur_K_init(end)*0.9;
Nocibur_periods(end+1) = Nocibur_periods(end);
Nocibur_K_init(end+1) = Nocibur_K_init(end)*0.9;
Nocibur_periods(end+1) = Nocibur_periods(end);
Nocibur_K_init(end+1) = Nocibur_K_init(end)*0.9;
Nocibur_periods(end+1) = Nocibur_periods(end);
Nocibur_K_init(end+1) = Nocibur_K_init(end)*0.9;
end
git_par.z_stop = sum(Nocibur_periods);
GIT_fun_2019;
GIT_fun_output;

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

tic;

% GIT_dir should be in MATLAB path (or in same directory as input file)
GITMATLAB_path= genpath('C:\Users\pietro\Desktop\GIT_2019');
addpath(GITMATLAB_path);
GIT_dir;
%% Parameter list to run GIT_fun 
git_par.prad0 = 1e9;    
git_par.parallel = 0;            % 1 parallel ; 0 single processor 
git_par.helical = 1;             % 0 for planar; 1 for helical
git_par.period_tapering = 0;    % 0 for constant period, 1 for varying period, -1 for predetermined period
git_par.verbose = 0;                 % plot phase space during run
git_par.outputfrequency = 1;   % # of steps between command window status report 
git_par.plot_phase=0;
git_par.fieldcontourplot = 0;      % plot central slice transverse profile contour plot at the end of each Genesis run
git_par.writevid=0;                     % write a video file with phasespace evolution

%% Different TESSA optimization schemes
git_par.lps_predictor = 0;       % predicts phase space to decide tapering step. only works for NWIG small - obsolete
git_par.correct_theta = 0;       % every drift section a phase shifter pushes the trapped phasespace in the center of the resonant bucket
git_par.correct_R56 =   1;       % Correct for artificial r56 introduced by Genesis in drifts 

git_par.psi_0 = 0.9;               % resonant phase
git_par.fac0 = 1;                   % For uniform distribution use fac0 = 1 psi_0 = 0.9 
                                            % For parabolic distribution use fac0 = 0.95 psi_0 = 1.3      

git_par.offset = 0;                      % slippage control parameter. should be -1.
 
git_par.quadrupole_lattice = 2;   % 1 Use quadrupole lattice from magin file
                                                   % 0 Keep a constant beam size using xkx and xky
                                                   % 2 Use long undulator sections
                               
git_par.rxbeam = 30.00e-6;   % Matched beam radius
git_par.rybeam = 30.00e-6;   % Matched beam radius

%% GENESIS parameters
git_par.itgaus = 2;     % 0 parabolic ; 1 gaussian ; 2 uniform transverse profiles
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
git_par.ncar = 121;
git_par.dgrid = 3e-3;
git_par.nscz = 0;
git_par.nscr = 0;
git_par.pulselength = 180e-6;    % Laser pulse length
git_par.beamfile = 1;
git_par.radfile = 1;

%% parameters used for period tapering 
git_par.gap = 0.003;
git_par.Brem = 1.5;

%% parameters for transverse matching. Calculated separately.   
%%%48T/m 227T/m parameter
% beamsize_pbmatch;

git_par.rxbeam_match=94.9e-06;
git_par.rybeam_match=115.23e-06;

git_par.alphax=   -0.822;
git_par.alphay=    1.577;
git_par.quadgradient = 0;
git_par.xkx = 0.5;
git_par.xky = 0.5;
git_par.lquad = 0;
git_par.fododrift = 0;

%% Run GIT_fun (should be in same directory or in Matlab path)
recirculated = 0;

%% Tapered undulator
%beamfile and radfile
git_par.nslice=800;         %number of slices in Genesis
git_par.curpeak = 1000;     %current peak
git_par.curlen = 90e-6;     % Electron beam pulse length
curpeakpos=git_par.nslice/2; %position of current peak
radshift=curpeakpos-120;    %position of radiation peak
beamtype=0; %beamtype==0: Rectangular, beamtype==1:Gaussian
write_radfile;
write_beamfile;
p_rad_0=p_rad; %for plotting 

%% prebuncher
jind = 1;
for iy = 1:1:30
    for ix = 1:1:20
        git_par.prad0 = iy*1e9;
        R56buncher = ix*1e-6;      
prebunched = 0;
git_par.nwig=53;
git_par.quadrupole_lattice=2;
git_par.z_stop =.032*53;          % Simulation length
git_par.period_tapering = 0;
git_par.psi_0 = 0;            
totalmagfilename=['mg_pb.mag.in'];
copyfile(totalmagfilename,[dirname,'mg.mag.in']);
git_par.awd=0;
readin_totalmagfile;
git_par.drift = 45;                  % drift length in undulator period units 
GIT_fun_2019;
p_rad_1=power(:,end);
phaseshift = 1;
prop_part;
%plot_part_phase;
    y(jind) =git_par.prad0; 
    x(jind) = R56buncher;
    outres(jind) = abs(bcomplex);
    jind = jind + 1;
    end
end
%%
[xq,yq] = meshgrid(1e-6:1e-6:18e-6,1e9:1e9:30e9);
vq = griddata(x,y,outres,xq,yq);
contourf(vq);
xlabel('R56 (microns)')
ylabel('Seed Power (GW)')
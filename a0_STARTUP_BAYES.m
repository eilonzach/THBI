%% script to put various things on path neded for inversion 

% turn off warning about name conflict with matlab builtin isstring
warning('off','MATLAB:dispatcher:nameConflict');
% set time zone:
setenv('TZ','America/Los_Angeles')

basedir = '/Users/zeilon/Documents/MATLAB/';
bayesdir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/';

% path of all inversion main functions
addpath(bayesdir)
% path of all inversion sub functions
addpath([bayesdir,'functions'])
% path to fast spline func.
addpath('/Users/zeilon/Dropbox/MATLAB/lib/fastBSpline'); 
% path to propagator matrix running dir.
addpath([bayesdir,'matlab_to_propmat']); 
% path to mineos running dir.
addpath([basedir,'matlab_to_mineos']); 
% path to seiz models
addpath('~/Dropbox/MATLAB/lib/seizmo/models')
% path to HV kernels functions
addpath('~/Work/codes/HV_Tanimoto/matlab_to_HVkernel')
% path to Rayleigh wave dispersion curve dir
addpath('/Users/zeilon/Dropbox/MATLAB/seis_tools/surface_waves'); 
% path to gaussfit dir
addpath('~/Dropbox/MATLAB/lib/gaussfit'); 

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict');

global prem_anisotropic prem_isotropic
prem_anisotropic = prem_perfect('SPVW',0.25);
prem_isotropic = prem;
%% script to put various things on path neded for inversion 

% turn off warning about name conflict with matlab builtin isstring
warning('off','MATLAB:dispatcher:nameConflict');

basedir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/';
% path of all inversion main functions
addpath(basedir)
% path of all inversion sub functions
addpath([basedir,'functions'])
% path to fast spline func.
addpath('/Users/zeilon/Dropbox/MATLAB/lib/fastBSpline'); 
% path to propagator matrix running dir.
addpath([basedir,'matlab_to_propmat']); 
% path to mineos running dir.
addpath([basedir,'matlab_to_mineos']); 
addpath('~/Documents/MATLAB/seizmo/models/')
% path to Rayleigh wave dispersion curve dir
addpath('/Users/zeilon/Dropbox/MATLAB/seis_tools/surface_waves'); 

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict');

global prem_anisotropic prem_isotropic
prem_anisotropic = prem_perfect('SPVW',0.25);
prem_isotropic = prem;
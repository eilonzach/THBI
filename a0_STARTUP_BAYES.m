%% script to put various things on path neded for inversion 

% turn off warning about name conflict with matlab builtin isstring
warning('off','MATLAB:dispatcher:nameConflict');
% set time zone:
setenv('TZ','America/Los_Angeles')

basedir = '/Users/zeilon/Documents/MATLAB/'; % change these to your values
bayesdir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/'; % change these to your values

% You'll likely need to get all of these functions - email me if not
% obvious. Some turn out not to be used any more, so comment out and then
% let me know if they come up and you need them. 

% path of all inversion main functions
addpath(bayesdir)
% path of all inversion sub functions
addpath([bayesdir,'functions']) % inside this folder
% path to fast spline func.
addpath('/Users/zeilon/Dropbox/MATLAB/lib/fastBSpline'); % google "fastBSspline MATLAB" for this
% path to propagator matrix running dir.
addpath([bayesdir,'matlab_to_propmat']); % provided
% path to mineos running dir.
addpath([basedir,'matlab_to_mineos']); % provided
% path to seiz models
addpath('~/Dropbox/MATLAB/lib/seizmo/models') % Google "MATLAB seizmo" for this
% path to HV kernels functions
addpath('~/Work/codes/HV_Tanimoto/matlab_to_HVkernel') % only if doing ellipticity ratios
% path to Rayleigh wave dispersion curve dir
addpath('/Users/zeilon/Dropbox/MATLAB/seis_tools/surface_waves'); % on my github
% path to gaussfit dir
addpath('~/Dropbox/MATLAB/lib/gaussfit'); % don't think used any more

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict');

global prem_anisotropic prem_isotropic
prem_anisotropic = prem_perfect('SPVW',0.25);
prem_isotropic = prem;
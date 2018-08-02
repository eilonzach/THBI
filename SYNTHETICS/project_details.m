%% Initialise
proj = struct('name','SYNTHETICS');
proj.dir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/SYNTHETICS';

%% directories 
% proj.rawdatadir = '/Volumes/data/THBI/US/STAsrawdat';
proj.STAinversions = proj.dir;
datN = '';

%% Add matguts to load data to the path
addpath('~/Dropbox/myPAPERS/2018_THBI_methods/codes/Figure_3/')

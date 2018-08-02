%% Initialise
proj = struct('name','LAB_tests');
proj.dir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/LAB_tests';

%% directories 
% proj.rawdatadir = '/Volumes/data/THBI/US/STAsrawdat';
proj.STAinversions = proj.dir;
datN = '';

%% Add matguts to load data to the path
addpath('~/Dropbox/myPAPERS/2018_THBI_methods/codes/Figure_3/')

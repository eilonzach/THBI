%% Initialise
proj = struct('name','US');
proj.dir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/US';

%% directories 
proj.rawdatadir = '/Volumes/data/THBI/US/STAsrawdat';
proj.STAinversions = '/Volumes/data/THBI/US/STAsinv';

%% Add matguts to load data to the path
addpath([proj.dir,'/matguts']);

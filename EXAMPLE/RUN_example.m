close all
clear all
%% Setup
run('../a0_STARTUP_BAYES.m')

proj = struct('name','EXAMPLE');
proj.dir = [fileparts(mfilename('fullpath')),'/'];
proj.STAinversions = [proj.dir,'inversion_results/'];
save([proj.dir,'project_details.mat'],'proj')

wd = pwd; addpath(wd);
cd(proj.dir);


%% specify details of this run
generation = 0; % generation of solution and data processing
gc = '';
BWclust = '';
STAMP = 'EXAMPLE';

onesta = '';

%% put parameters in place 
global run_params

% you can obviously adapt these (likely in some loop over stations etc.) 
% to be appropriate for your dataset
run_params.projname = proj.name; % from above
run_params.sta = 'test_station'; % name of station
run_params.nwk = 'test_nwk'; % name of network
run_params.gc = gc; % great circle distance of body wave data, if relevant
run_params.BWclust = BWclust; % cluster of BW data, if relevant
run_params.datN = generation; % processing iteration, if relevant
run_params.STAMP = STAMP; % NEED - some identifier for this inversion run
run_params.overwrite = 1; % do you want to overwrite previous results?

%% Run it
MASTER_par;

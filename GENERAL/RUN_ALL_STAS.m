close all
clear all
%% Setup
proj = struct('name','NWUS');
proj.dir = ['~/Documents/MATLAB/BayesianJointInv/',proj.name];
wd = pwd;
addpath(wd);
cd(proj.dir);

%% load project, station, and request details and request details
load([proj.dir,'/project_details.mat']);
load([proj.infodir,'stations.mat']);

%% specify details of this run
generation = 30; % generation of solution and data processing
gc = 'all';

STAMP = 'AGU18v1';

overwrite = true;

%% put parameters in place for running all stations
global run_params

run_params.projname = proj.name;
run_params.gc = gc;
run_params.datN = generation;
run_params.STAMP = STAMP;
run_params.overwrite = overwrite;

%% ==================  LOOP OVER STATIONS IN DB  ================== 
for is = 10:stainfo.nstas
    fprintf('\n\n'); for i=1:3, for j=1:40, fprintf('**'); end; fprintf('*\n'); end; fprintf('\n');
    
    fprintf('STATION: %s\n',stainfo.stas{is})
    fprintf('NETWORK: %s\n\n',stainfo.nwk{is})
    run_params.sta = stainfo.stas{is};
    run_params.nwk = stainfo.nwk{is};
    
    execute_MASTER_par
end

function execute_MASTER_par
    global run_params
    try
    MASTER_par;
    end
end
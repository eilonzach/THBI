close all
clear all
%% Setup
run('/Users/zeilon/Documents/MATLAB/BayesianJointInv/a0_STARTUP_BAYES.m')

proj = struct('name','NWUS');
proj.dir = ['~/Documents/MATLAB/BayesianJointInv/',proj.name];
wd = pwd; addpath(wd);
cd(proj.dir);

%% load project, station, and request details and request details
try
    load([proj.dir,'/project_details.mat']);
catch
    run([proj.dir,'/project_details.m']);
end
load([proj.infodir,'stations.mat']);

%% specify details of this run
generation = 30; % generation of solution and data processing
gc = 1;
BWclust = 1;
STAMP = 'NWUS_CCP';

onesta = '';

overwrite = true;

%% put parameters in place for running all stations
global run_params

run_params.projname = proj.name;
run_params.gc = gc;
run_params.BWclust = BWclust;
run_params.datN = generation;
run_params.STAMP = STAMP;
run_params.overwrite = overwrite;

%% ==================  LOOP OVER STATIONS IN DB  ================== 
for is = 100:stainfo.nstas % got to BMO, UO = is 89
    
    if exist('onesta') && ~isempty(onesta)
        if ~strcmp(stainfo.stas{is},onesta), continue; end
    end
    
    fprintf('\n'); for i=1:3, for j=1:40, fprintf('**'); end; fprintf('*\n'); end; fprintf('\n');
    
    fprintf('STATION: %s\n',stainfo.stas{is})
    fprintf('NETWORK: %s\n\n',stainfo.nwk{is})
    run_params.sta = stainfo.stas{is};
    run_params.nwk = stainfo.nwk{is};
    
    % do the work (and make all the Mineos files) in a workdir
    if exist('workdir','dir')~=7, mkdir('workdir'); end
%     cd('workdir')
    cd(proj.dir)
    
    execute_MASTER_par

    cd(proj.dir)
end

function execute_MASTER_par
    global run_params
    try
    MASTER_par;
    catch
%         1;
    end
end
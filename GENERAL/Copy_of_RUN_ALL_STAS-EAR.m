close all
clear all
%% Setup
proj = struct('name','EARdb');
proj.dir = ['~/Work/EastAfrica/',proj.name];
wd = pwd; addpath(wd);
cd(proj.dir);

%% load project, station, and request details and request details
load([proj.dir,'/project_details.mat']);
load([proj.infodir,'stations.mat']);

%% specify details of this run
generation = 30; % generation of solution and data processing
gc = 1;
BWclust = 1;
STAMP = 'TEI19_BWSW';

onesta = '';

overwrite = false;

%% put parameters in place for running all stations
global run_params

run_params.projname = proj.name;
run_params.gc = gc;
run_params.BWclust = BWclust;
run_params.datN = generation;
run_params.STAMP = STAMP;
run_params.overwrite = overwrite;

%% ==================  LOOP OVER STATIONS IN DB  ================== 
for is = 59:stainfo.nstas % got to BMO, UO = is 89
    
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
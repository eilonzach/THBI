%% Script to download the data for all stations for this project
clear all

proj = struct('name','AFR');
proj.dir = ['~/Documents/MATLAB/BayesianJointInv/',proj.name];
wd = pwd;
addpath(wd);

%% load project, station, and request details and request details
load([proj.dir,'/project_details']);
load([proj.infodir,'stations']);
load([proj.infodir,'/data_request_details']);

return
%% Request data    
if exist([proj.rawdatadir,'request_files/'],'dir')~=7 
    mkdir([proj.rawdatadir,'request_files/']); 
end

for is = 1:stainfo.nstas
    [reqfile] = evdata_WAVEFORMS_breqfast(stainfo.stas{is},stainfo.nwk{is},true,false,request_details_all);
    for ii = 1:length(reqfile)
        movefile(reqfile{ii},[proj.rawdatadir,'request_files/']);
    end
end

return
%% Download and process data
for is = 1:stainfo.nstas
    [~, datafile] = evdata_WAVEFORMS_breqfast(stainfo.stas{is},stainfo.nwk{is},false,true,request_details_all);
	movefile([datafile,'.mat'],proj.rawdatadir);
end
return

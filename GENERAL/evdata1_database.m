%% Script to establish a database of stations and events for body wave study
clear all
ifrunall = true;

proj = struct('name','NWUS');
proj.dir = ['~/Documents/MATLAB/BayesianJointInv/',proj.name];

%% Station parameters
sta_latlims = [40 55]; % [min_lat max_lat] for stations
sta_lonlims = [-130 -85]; % [min_lon max_lon] for stations
sta_chans = 'BH*,HH*'; % channel codes to search for
starttime = '1970-01-01 00:00:00';
startbytime = '2010-01-01 00:00:00';
min_longevity_yrs = 10;

%% Event parameters
mag_lims = [5.7 7.4];
dep_lims = [0 1000]; % set to [0 1000] by default (km)
gc_lims  = [30 75];
% startafter = '1990-01-01 00:00:00'; % earliest evtime (yyyy-mm-dd HH:MM:SS)

%% Data parameters
phases = {'P','S'};
samprate = 40;
datwind =  [-100 100];

%% ID for IRIS DMC request
IRIS_ID = 'zeilon';


%% GET TO WORK
wd = pwd;
addpath('matguts');

%% Make directory structure
% main database directory
proj.dir = regexprep(proj.dir,'~',getenv('HOME'));
if ~strcmp(proj.dir(end),'/'),proj.dir = [proj.dir,'/']; end
if exist(proj.dir,'dir')~=7, mkdir(proj.dir); end
% info files directory 
proj.infodir = [proj.dir,'INFO/'];
if exist(proj.infodir,'dir')~=7, mkdir(proj.infodir); end
% response files directory 
proj.respdir = [proj.dir,'INFO/RESP/'];
if exist(proj.respdir,'dir')~=7, mkdir(proj.respdir); end
% data files directory 
proj.rawdatadir = ['/Volumes/data/THBI/',proj.name,'/STAsrawdat/'];
proj.STAinversions = ['/Volumes/data/THBI/',proj.name,'/STASinv/'];
if exist(proj.rawdatadir,'dir')~=7, mkdir(proj.rawdatadir); end
if exist(proj.STAinversions,'dir')~=7, mkdir(proj.STAinversions); end

% save project details
save([proj.dir,'project_details.mat'],'proj');

%% Add matguts to load data to the path
addpath([proj.dir,'/matguts']);


cd(proj.dir)

if ~ifrunall
    return 
end

%% Write request details information
request_details_all = struct(...
        'phases',	{phases},... % need double {{ }}
        'gclims',   gc_lims,...
        'maglims',  mag_lims,...
        'samprate', samprate,...
        'datwind',  datwind);

save([proj.infodir,'/data_request_details.mat'],'request_details_all');

if ~ifrunall
    return
end

%% Get station + channel information
javaaddpath('/Users/zeilon/Documents/MATLAB/IRIS-WS-2.0.15.jar')
% save station request info
stations_request = struct('lat_lims',sta_latlims,'lon_lims',sta_lonlims,'chans',sta_chans,...
                          'starttime',starttime,'startbytime',startbytime,'min_longevity_yrs',min_longevity_yrs);

% grab stations
stations_IRIS = irisFetch.Stations('station','*','*','*','BH?',...
    'boxcoordinates',[sta_latlims,sta_lonlims],'StartAfter',starttime,'StartBefore',startbytime);
%only include stations satisfying longevity
starter = datenum({stations_IRIS.StartDate}');
ender = datenum({stations_IRIS.EndDate}'); ender(ender>now)= now;
stations_IRIS = stations_IRIS((ender - starter)/365.25 >= min_longevity_yrs)


stainfo = struct('stas',{{stations_IRIS.StationCode}'},...
                 'nwk',{{stations_IRIS.NetworkCode}'},...
                 'slats',[stations_IRIS.Latitude]',...
                 'slons',[stations_IRIS.Longitude]',...
                 'selevs',[stations_IRIS.Elevation]',...
                 'ondate',datenum({stations_IRIS.StartDate}'),...
                 'offdate',datenum({stations_IRIS.EndDate}'),...
                 'ondate_str',{{stations_IRIS.StartDate}'},...
                 'offdate_str',{{stations_IRIS.EndDate}'},...
                 'nstas',length(stations_IRIS));  
             
[stainfo] = stainfo_unique(stainfo);

% parse channels             
chans = cell(stainfo.nstas,3);
chandips = nan(stainfo.nstas,3);
chanazs = nan(stainfo.nstas,3);
nchans = zeros(stainfo.nstas,1);
for is = 1:stainfo.nstas
    nchan = length(stations_IRIS(is).Channels);
    tempchans = cell(1,nchan);
    tempdips = nan(1,nchan);
    tempazs = nan(1,nchan);
    for ic = 1:nchan
        tempchans(ic) = {stations_IRIS(is).Channels(ic).ChannelCode};
        tempdips(ic) = stations_IRIS(is).Channels(ic).Dip;
        tempazs(ic) = stations_IRIS(is).Channels(ic).Azimuth;
    end
    [stachans,indch] = unique(tempchans);
    nchans(is) = length(stachans);
    chans(is,1:nchans(is)) = stachans;
    chandips(is,1:nchans(is)) = tempdips(indch);
    chanazs(is,1:nchans(is)) = tempazs(indch);
end
chandips(cellfun('isempty',chans)) = nan;
chanazs(cellfun('isempty',chans)) = nan;
stainfo.nchans = nchans;
stainfo.chans = chans;
stainfo.chandips = chandips;
stainfo.chanazs = chanazs;

                        
save([proj.infodir,'/stations'],'stainfo','stations_IRIS','stations_request');

if ~ifrunall
    return
end

%% Response SAC_PZ files
% Build and send BREQFAST request file for dataless seed                     
addpath('~/Dropbox/MATLAB/seis_tools/breqfasting/');
breq_fast_request([proj.name,'_dataless'],IRIS_ID,{stations_IRIS.StationCode}','*',{stations_IRIS.NetworkCode}','',{stations_IRIS.StartDate}',{stations_IRIS.EndDate}','dataless_SEED',[proj.name,'_dataless_request']);
movefile([proj.name,'_dataless_request'],proj.respdir);
return 
%% ======= WAIT FOR NOTIFICATION THAT DATALESS SEED IS ON SERVER ======= %%
% Download and process dataless seed  
breq_fast_dataless_PZprocess([proj.name,'_dataless'],IRIS_ID,proj.respdir,{'BH*','HH*'},1 )     
return
                        
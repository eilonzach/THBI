clear all
close all

datadir = '/Volumes/zeilon/USdata';

if exist(datadir,'dir')~=7, mkdir(datadir); end

request_details_all = struct(...
        'phases',	{{'P','S'}},... % need double {{ }}
        'gclims',   [30 75],...
        'maglims',  [5.8 7.3],...
        'samprate', 40,...
        'datwind',  [-150 150]);


evdata1_stations
stainfo_master = stainfo;

% request data
for is = 1:length(stainfo_master)
    dfile = dir([datadir,'/dat_',stainfo_master(is).StationCode,'_*_30to75.mat']);
    if ~isempty(dfile), continue, end
    fprintf('\n============================================\n')
    fprintf('Requesting %s %s...',stainfo_master(is).StationCode,stainfo_master(is).NetworkCode)
    reqfile = evdata2_WAVEFORMS_breqfast(stainfo_master(is).StationCode,...
                                         stainfo_master(is).NetworkCode,...
                                         true,false,request_details_all);
    for ip = 1:length(reqfile)
        movefile(reqfile{ip},datadir)
    end
end


pause(3*60*60)
% process data
for is = 1:length(stainfo_master)
    dfile = dir([datadir,'/dat_',stainfo_master(is).StationCode,'_*_30to75.mat']);
    if ~isempty(dfile), continue, end
    fprintf('\n============================================\n')
    fprintf('Downloading %s %s...',stainfo_master(is).StationCode,stainfo_master(is).NetworkCode)
    try 
    [~,datafile] = evdata2_WAVEFORMS_breqfast(stainfo_master(is).StationCode,...
                                         stainfo_master(is).NetworkCode,...
                                         false,true,request_details_all);
    movefile([datafile,'.mat'],datadir)                                     
    catch me
        fprintf('... NO DATA?!\n')
    end
end

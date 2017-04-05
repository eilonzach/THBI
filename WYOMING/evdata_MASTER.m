clear all
close all

if exist('ALLrequests','dir')~=7, mkdir('ALLrequests'); end


evdata1_stations
stainfo_master = stainfo;

% request data
for is = 1:length(stainfo_master)
    if ~strcmp(stainfo(is).StationCode,'E19A'), continue; end
    fprintf('\n============================================\n')
    fprintf('Requesting %s %s...',stainfo_master(is).StationCode,stainfo_master(is).NetworkCode)
    reqfile = evdata2_WAVEFORMS_breqfast(stainfo_master(is).StationCode,...
                                         stainfo_master(is).NetworkCode,...
                                         true,false);
    for ip = 1:length(reqfile)
        movefile(reqfile{ip},'ALLrequests')
    end
end
return

pause(3*60*60)
% process data
for is = 1:length(stainfo_master)
    fprintf('\n============================================\n')
    fprintf('Downloading %s %s...',stainfo_master(is).StationCode,stainfo_master(is).NetworkCode)
    try 
    [~,datafile] = evdata2_WAVEFORMS_breqfast(stainfo_master(is).StationCode,...
                                         stainfo_master(is).NetworkCode,...
                                         false,true);
    movefile(datafile,'DATA')                                     
    catch
        fprintf('... NO DATA?!\n')
    end
end

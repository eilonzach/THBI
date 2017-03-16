%% Script to grab station information for Wyoming craton project

lonlim = [-113 -93];
latlim = [40 50];

starttime = '1970-01-01 00:00:00';
startbytime = '2005-01-01 00:00:00';
min_longevity_yrs = 10;

%% Get to work
javaaddpath('/Users/zeilon/Documents/MATLAB/IRIS-WS-2.0.15.jar')
addpath('plotting')

% grab stations
stainfo = irisFetch.Stations('station','*','*','*','BH?',...
    'boxcoordinates',[latlim,lonlim],'StartAfter',starttime,'StartBefore',startbytime);
%only include stations satisfying longevity
stainfo = stainfo((datenum({stainfo.EndDate}') - datenum({stainfo.StartDate}'))/365.25 >= min_longevity_yrs)

figure(1), clf
plot([stainfo.Longitude],[stainfo.Latitude],'o')
text([stainfo.Longitude]+0.1,[stainfo.Latitude]+0.2,{stainfo.StationName})
text([stainfo.Longitude]+0.1,[stainfo.Latitude],{stainfo.StationCode})
text([stainfo.Longitude]+0.1,[stainfo.Latitude]-0.2,{stainfo.NetworkCode})
add_state_boundaries(gca,latlim,lonlim)
set(gca,'ylim',latlim,'xlim',lonlim)


% Decent stations:
% 'LKWY, US' -Yellowstone Lake, Wyoming, USA
% 'RWWY, IW' -  Rawlins, Wyoming, USA
% 'RSSD, IU' - South Dakota
% 'BW06, US' - 	Boulder Array Site 6 (Pinedale Array Site 6), Wyoming, USA
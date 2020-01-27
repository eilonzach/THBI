function [ periods, phVs ] = Lph_dispcurve_latlon(ilat,ilon,transT)
% [ periods, phVs ] = Lph_dispcurve_latlon(ilat,ilon,transT)
% 
%  Function to obtain an ambient noise Love wave dispersion curve at any
%  lat/lon point . Each dispersion curve is interpolated from individual
%  phase velocity maps at the frequencies of interest. If the requested
%  lat/lon point is outside the grid of data, this function will return
%  phV_freq = nan.

if nargin<3 || isempty(transT)
    transT = 33;
end

seismoddir = '/Volumes/data/models_seismic/';
if ~exist(seismoddir,'dir')
    try
        seismoddir = '/Volumes/eilon_data/models_seismic/';
    catch
        error('NO SEISMOD DIR FOUND');
    end
end

%% EQ data:

% none, for now

%% AN data:
% EKSTROM
datadir = [seismoddir,'US_LOVE_ANT_phV_EKSTROM/'];
ANperiods = [6:2:12,15:5:40]';
ANphV_period = disp_curve_AN_latlon(ANperiods,ilat,ilon,datadir);

% %% composite
% p1 = ANperiods(ANperiods<min(Eperiods));
% v1 = ANphV_period(ANperiods<min(Eperiods));
% 
% p3 = Eperiods(Eperiods>transT);
% v3 = EphV_period(Eperiods>transT);
% 
% p2 = [min(Eperiods):2:transT]';
% v2 = mean([interp1(Eperiods,EphV_period,p2),interp1(ANperiods,ANphV_period,p2)],2);
% 
% periods = [p1;p2;p3];
% phVs    = [v1;v2;v3];
periods = ANperiods;
phVs = ANphV_period;

%% delete any nans
kill = isnan(phVs);
periods(kill) = [];
phVs(kill) = [];


end


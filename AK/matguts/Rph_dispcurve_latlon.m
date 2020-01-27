function [ periods, phVs ] = Rph_dispcurve_latlon(ilat,ilon,transT)
% [ periods, phVs ] = Rph_dispcurve_latlon(ilat,ilon,transT)
% 
%  Function to obtain a composite dispersion curve at any lat/lon point by
%  combination of Earthquake and Ambient noise dispersion curves. Each
%  dispersion curve is interpolated from individual phase velocity maps at
%  the frequencies of interest. If the requested lat/lon point is outside
%  the grid of data, this function will return phV_freq = nan.

if nargin<3 || isempty(transT)
    transT = 33;
end

%% EQ data:

ddir = '/Volumes/data/models_seismic/US_RAYLEIGH_EQ_phV_DALTON/';
Eperiods = [25,40,50,60,80,100,120,140,180]';
EphV_period = disp_curve_EQ_latlon(Eperiods,ilat,ilon,ddir);
[Eperiods,iE] = sort(Eperiods);
EphV_period = EphV_period(iE);

%% AN data:
% try Shen and Ritzwoller, 2016
datadir = '/Volumes/data/models_seismic/US_RAYLEIGH_ANT_phV_SHEN/';
ANperiods = [8:2:32,36,40]';
ANphV_period = disp_curve_AN_latlon(ANperiods,ilat,ilon,datadir);

%% composite
p1 = ANperiods(ANperiods<min(Eperiods));
v1 = ANphV_period(ANperiods<min(Eperiods));

p3 = Eperiods(Eperiods>transT);
v3 = EphV_period(Eperiods>transT);

p2 = [min(Eperiods):2:transT]';
v2 = mean([interp1(Eperiods,EphV_period,p2),interp1(ANperiods,ANphV_period,p2)],2);

periods = [p1;p2;p3];
phVs    = [v1;v2;v3];

%% delete any nans
kill = isnan(phVs);
periods(kill) = [];
phVs(kill) = [];


end


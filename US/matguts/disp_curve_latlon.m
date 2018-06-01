function [ periods, phVs ] = disp_curve_latlon(ilat,ilon,transT)
% [ periods, phVs ] = disp_curve_latlon(ilat,ilon)
% 
%  Function to obtain a composite dispersion curve at any lat/lon point by
%  combination of Earthquake and Ambient noise dispersion curves. Each
%  dispersion curve is interpolated from individual phase velocity maps at
%  the frequencies of interest. If the requested lat/lon point is outside
%  the grid of data, this function will return phV_freq = nan.

if nargin<3 || isempty(transT)
    transT = 33;
end

[ freqs, Eperiods ] = get_freqs;
EphV_period = disp_curve_EQ_latlon(freqs,ilat,ilon);
[Eperiods,iE] = sort(Eperiods);
EphV_period = EphV_period(iE);

%% AN data:
ANperiods = [8:2:32,36,40]';
ANphV_period = disp_curve_AN_latlon(ANperiods,ilat,ilon);

%% composite
p1 = ANperiods(ANperiods<min(Eperiods));
v1 = ANphV_period(ANperiods<min(Eperiods));

p3 = Eperiods(Eperiods>transT);
v3 = EphV_period(Eperiods>transT);

p2 = [min(Eperiods):2:transT]';
v2 = mean([interp1(Eperiods,EphV_period,p2),interp1(ANperiods,ANphV_period,p2)],2);

periods = [p1;p2;p3];
phVs    = [v1;v2;v3];


end


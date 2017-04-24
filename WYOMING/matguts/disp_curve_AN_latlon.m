function [ phV_period ] = disp_curve_AN_latlon( periods,ilat,ilon,datadir )
% [ phV_period ] = disp_curve_AN_latlon( freqs,ilat,ilon )
%  Function to obtain a dispersion curve at any lat/lon point by
%  interpolation of individual ambient noise phase velocity maps at the
%  periods of interest. If the requested lat/lon point is outside the grid
%  of data, this function will return phV_freq = nan.

if nargin < 4 || isempty(datadir)
    datadir = '~/Work/data/models_seismic/US_RAYLEIGH_PHASE_ANT_VEL_SHEN/'; % need final slash
end

phV_period = nan(length(periods),1);

for ii = 1:length(periods)
    [lats,lons,phVs] = load_ANphV_data( periods(ii),datadir );
    phV_period(ii) = griddata(lons,lats,phVs,ilon,ilat);
end




end


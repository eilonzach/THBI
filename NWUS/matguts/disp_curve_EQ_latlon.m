function [ phV_freq ] = disp_curve_EQ_latlon( freqs,ilat,ilon,datadir )
% [ phV_freq ] = disp_curve_latlon( freqs,ilat,ilon )
%  Function to obtain a dispersion curve at any lat/lon point by
%  interpolation of individual phase velocity maps at the frequencies of
%  interest. If the requested lat/lon point is outside the grid of data,
%  this function will return phV_freq = nan.

if nargin < 4 || isempty(datadir)
	datadir = '/Volumes/data/models_seismic/US_RAYLEIGH_EQ_phV_DALTON/'; % need final slash
end

phV_freq = nan(length(freqs),1);

for ii = 1:length(freqs)
    [lats,lons,phVs] = load_phV_data( freqs(ii),datadir );
    phV_freq(ii) = griddata(lons,lats,phVs,ilon,ilat);
end




end


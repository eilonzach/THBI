function [ err ] = error_curve_latlon( freqs,ilat,ilon,errdir )
% [ err ] = error_curve_latlon( freqs,ilat,ilon,errdir )
%  Function to obtain the errors in phase velocity at any lat/lon point by
%  interpolation of individual phase velocity error maps at the frequencies
%  of interest. If the requested lat/lon point is outside the grid of data,
%  this function will return err = nan.

if nargin < 4 || isempty(errdir)
	errdir = '~/Dropbox/Dave_Li_phV/Errors/'; % need final slash
end

err = nan(length(freqs),1);

for ii = 1:length(freqs)
    [lats,lons,errs] = load_phV_error( freqs(ii),errdir );
    err(ii) = griddata(lons,lats,errs,ilon,ilat);
end




end


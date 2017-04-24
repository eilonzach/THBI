function [ lat,lon,phV, latgrid, longrid,phVgrid ] = load_ANphV_data( period,datadir )
% [ lat,lon,phV ] = load_ANphV_data( period,datadir )
%   Function to load a ambient noise phase velocity map, given the
%   frequency of interest. Assumes that file names are in the format
%   "PERIOD_ANT.vel", with coluns of lon, lat, phV.

if nargin<2 || isempty(datadir)
    datadir = '~/Work/data/models_seismic/US_RAYLEIGH_PHASE_ANT_VEL_SHEN/'; % need final slash
end

% name the file
datafile = [num2str(period),'_ANT.vel'];

% load the data
dat = load([datadir,datafile]);

% parse the data into columns
lon = dat(:,1); lon = mod(lon+180,360)-180; % make -180 to 180
lat = dat(:,2); 
phV = dat(:,3);


if nargout>3;
% grid the data
unqlat = unique(lat);
unqlon = unique(lon);
[longrid,latgrid] = meshgrid(unqlon,unqlat);
phVgrid = griddata(lon,lat,phV,longrid,latgrid);
end

end


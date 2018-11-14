function [ lat,lon,phV, latgrid, longrid,phVgrid ] = load_phV_data( freq,datadir )
% [ lat,lon,phV ] = load_phV_data( freq,datadir )
%   Function to load a phase velocity map, given the frequency of interest.
%   Assumes that files are in the format "gridvalue.velFFFkern.dat", where
%   "FFF" is a three-character value for the frequency (e.g. "006") in mHz.
%   Input frequency should be in Hz.

if nargin<2 || isempty(datadir)
	datadir = '/Volumes/data/models_seismic/US_RAYLEIGH_EQ_phV_DALTON/'; % need final slash
end

% name the file
phVfile = sprintf('gridvalue.vel%03.0fkern.dat',freq*1000);

% load the data
fid = fopen([datadir,phVfile],'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

% parse the data into columns
lon = A{1};
lat = A{2};
phV = A{3};

if nargout>3;
% grid the data
unqlat = unique(lat);
unqlon = unique(lon);
[longrid,latgrid] = meshgrid(unqlon,unqlat);
phVgrid = griddata(lon,lat,phV,longrid,latgrid);
end

end


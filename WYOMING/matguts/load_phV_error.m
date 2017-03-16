function [ lat,lon,err, latgrid, longrid,errgrid ] = load_phV_error( freq,errordir )
% [ lat,lon,err ] = load_phV_data( freq,datadir )
%   Function to load a map of error in phase velocity, given the frequency
%   of interest. Assumes that files are in the format
%   "gridvalue.velFFFkern-covar.dat", where "FFF" is a three-character
%   value for the frequency (e.g. "006") in mHz. Input frequency should be
%   in Hz.

if nargin<2 || isempty(errordir)
    errordir = '~/Dropbox/Dave_Li_phV/Errors/'; % need final slash
end

% name the file
errfile = sprintf('gridvalue.vel%03.0fkern-covar.dat',freq*1000);

% load the data
fid = fopen([errordir,errfile],'r');
A = textscan(fid,'%f %f %f');
fclose(fid);

% parse the data into columns
lon = A{1};
lat = A{2};
err = A{3};

if nargout>3;
% grid the data
unqlat = unique(lat);
unqlon = unique(lon);
[longrid,latgrid] = meshgrid(unqlon,unqlat);
errgrid = griddata(lon,lat,err,longrid,latgrid);
end

end


function writeRAYSUMgeomfile( bazs,rayps,geomfile )
%writeRAYSUMgeomfile( parms,geomfile )
%   
% Function to write geometry file for RAYSUM code
% 
% INPUTS:
%  bazs     - vector of back azimuths (deg)
%  rayps    - vector of ray parameters (s/m)
%  geomfile - file name for output file

N = length(bazs); % # of rays

rayps_ex = floor(log10(rayps)); % base 10 exponent of ray parameter
rayps_ma = rayps.*10.^(-rayps_ex); % mantissa of ray parameter


% write
fid = fopen(geomfile,'w');
fprintf(fid,'# Columns: back-azimuth (deg), slowness (s/m), N-S shift (meters\n');
fprintf(fid,'# north), E-W shift (meters east).\n');
for ii = 1:N
    fprintf(fid,'%5.1f %1.3fE%.0f 0. 0.\n',bazs(ii),rayps_ma(ii),rayps_ex(ii)) ;   
end
fclose(fid);







end


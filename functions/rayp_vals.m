function [ unique_rayps,irayps ] = rayp_vals( rayps )
%[ unique_rayps,irayps ] = rayp_vals( rayps )
% 
% Function to parse a vector of ray parameters and determine the clusters
% of ray parameters that are similar and those that are different. The code
% will return the average ray parameters in each cluster (unique_rayps) and
% also the indices from the original vector that point to these individual
% % clusters (irayps). N.B. clusters within 1.0 s/deg linkages.
% 
% Assumes rayp in s/deg, but will try to adjust if in s/km


if length(rayps)==1
    unique_rayps = rayps;
    irayps = 1;
    return
end
    
Z = linkage(rayps(:));
irayps = cluster(Z,'cutoff',1.0,'criterion','distance');
N = max(irayps);

unique_rayps = zeros(N,1);
for ii = 1:N
    unique_rayps(ii) = mean(rayps(irayps==ii));
end


end


function [ wts ] = spline_weight( splines,spzz,zkni )
% [ wts ] = spline_weight( splines,zz )
%   Function to consider a set of splines with weights defined at a set of
%   depths and then pick out the spline(s) that have more than 40%
%   sensitivity at that depth (this will be 1/2 splines). This is so that
%   we can decide which splines to edit after we add a node.

w_z = interp1(spzz,splines,zkni);
wts = find(w_z>.4);
if isempty(wts)
    wts = mindex(-w_z);
end


end


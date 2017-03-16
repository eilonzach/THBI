function [ spbasis,spwts,knots ] = make_splines_fraction( model,par,fknots,m_or_c,zi,vi )
%  [ spbasis,spwts,knots ] = make_splines_fraction( model,par,fknots,m_or_c,zi,vi )
%   
% Function to make splines within a current model at fractions within the
% mantle/crust layer. Fractions will be rounded to the nearest 1/1000.We
% will find weights for these splines by interpolating the basis onto
% vectors zi and vi. 
% 
% N.B.  to move spline knots but not change velocities, just ignore the new
% weightings, use new spbasis and knots positions. 


if strcmp(m_or_c,'crust')    
    minz = model.sedmparm.h;
    maxz = minz + model.crustmparm.h;
elseif strcmp(m_or_c,'mantle')
    minz = model.sedmparm.h + model.crustmparm.h;
    maxz = par.mod.maxz + model.selev;
end

fknots = round_level(fknots,0.001); fknots = fknots(:);
knots = minz + fknots*(maxz-minz);
allzknots = [repmat(minz,3,1);knots;repmat(maxz,3,1)];

sp = fastBSpline.lsqspline(allzknots,2,zi,vi); % interpolate onto current model

zz = unique([minz:par.mod.dz:maxz,maxz])';
spbasis = sp.getBasis(zz); 
spbasis = spbasis(:,2:end-1);                

spwts = sp.weights(2:end-1); % pull back out spline coeff's from the interpolation        



end


function rho  = mantle_vs2rho( Vs,Zkm )
% rho  = mantle_vs2rho( Vs,Zkm )
%   Empirical scaling of Vs to rho for mantle rocks. density is scaled from
%   Vs using empirical scaling over a range of upper mantle P,T, where the
%   rho and Vs values are computed for "basal Lherzolite" using the
%   calculator of Hacker and Abers 2016 at conditions betwen 50 and 300 km,
%   with a mantle potential temperature of 1300 C, geotherm of 18C/GPa and
%   P calculated as Z(km)/32. From these values, we compute a best fitting
%   scaling function for Vs/rho as a function of pressure along the
%   geotherm and show that for temperature heterogeneity of ±60 C the error
%   in computed rho is less than 0.3%. See empirical_VsRho.m

tf = 1.337 + ((175-Zkm)/125)*0.0141;

rho = Vs./tf;

end


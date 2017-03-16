function Vp  = sed_vs2vp( Vs )
% vp  = sed_vs2vp( vs )
%   Empirical scaling of Vs to Vp for sedimentary rocks. 
%   Equations from Shen and Ritzwoller (JGR, 2016) equation (1), based
%   on results of Brocher (BSSA 2005).

Vp  =  0.941 + 2.095*Vs - 0.821*Vs.^2 + 0.268*Vs.^3 - 0.0251*Vs.^4;

end
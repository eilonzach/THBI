function [ q ] = q_v_given_V( vp,V,theta )
%[ q ] = q_v_given_V( vp,V,theta )
%   probability of a change in velocity from V to vp (v-prime) based on
%   gaussian probablilty function with standard deviation theta
% 
%  See Bodin et al., 2016, equation D4 (appendix D)

q = (1./(theta*sqrt(2*pi))) * exp( -0.5*(vp - V).^2 ./(theta.^2));


end


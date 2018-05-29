function [ q ] = q_v_given_V( vp,V,theta )
%[ q ] = q_v_given_V( vp,V,theta )
%   probability of a change in velocity from V to vp (v-prime) based on
%   gaussian probablilty function with standard deviation theta
% 
%  See Bodin et al., 2016, equation D4 (appendix D)

% q = gaussian(vp-V,0,theta)./gaussian(0,0,theta); % compare to prob of staying put

q = exp(((vp-V)/theta)^2/2)/theta/(sqrt(2*pi));

end


function [ misfit2,e ] = xconv_misfit( Vobs,Hobs,Vpre,Hpre )
%[ misfit2,e ] = xconv_misfit( Vobs,Hobs,Vpre,Hpre )
% 
%  Function to calculate the cross-convolution misfit between observed
%  (obs) and predicted (pre) time series on two components, V, and H. This
%  misfit function was (first?) proposed by Menke and Levin, 2003 for SKS
%  inversions, and has more recently been applied to receiver functions
%  (cf. Bodin 2013) where V is the vertical response and H is the
%  horizontal (radial) response for the P-SV system. 
%  N.B. here we DO NOT!! use conv(A,B,'same') to return only the central 
%  part of the cross-conv that
%  is the same size as the data vector. From what I can tell of the maths
%  in e.g. Bodin 2016 Appendix B (e.g. the size of Mv,Mh) so does Thomas.
% 
% INPUTS:
%   Vobs - observed vertical (or component 1)
%   Hobs - observed radial (or component 2)
%   Vpre - predicted vertical (or component 1)
%   Hpre - predicted radial (or component 2)
%  All inputs should be column vectors
% 
% OUTPUT:
%   misfit2 - squared cross-convolution misfit (a scalar) such that
%              e =  conv(Vobs,Hpre) - conv(Hobs,Vpre)
%              misfit2 = e'*e;
%   e       - vector of misfit btwn convolved pairs (same length as data)
% 
% Z. Eilon 08/2016

% idiot proofing
Vobs = Vobs(:);
Hobs = Hobs(:);
Vpre = Vpre(:);
Hpre = Hpre(:);

e = conv(Vobs,Hpre,'full') - conv(Hobs,Vpre,'full');
% e = conv(Vobs,Hpre,'same') - conv(Hobs,Vpre,'same');

misfit2 = e'*e;

end


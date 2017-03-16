function [ invCe,detCe,Ce ] = calc_Ce( sigma,r,n )
%[ invCe,detCe,Ce ] = calc_Ce( sigma,r,n )
%   Function to calculate the inverse and determinant of the data error
%   correlation matrix, from Bodin et al., GJR, 2012 - appendix D1
%  N.B. we assume the "first type" of noise parametrisation
% 
% sigma and r are values that parametrise the data
% n is the number of data points

%% calc detCe
detCe = (sigma.^(2*n))*((1 - r^2).^(n-1));

%% calc invCe
r1 = zeros(1,n); 
r1(1) = 1 + r^2; 
r1(2) = -r;

c1 = zeros(n,1);
c1(1) = 1 + r^2; 
c1(2) = -r;

M = toeplitz(c1,r1);
M(1,1) = 1;
M(n,n) = 1;

invCe = M./(sigma^2  * (1 - r^2));

invCe = sparse(invCe);

%% calc Ce0
if nargin>2
    r1 = zeros(n,1);
    for ii = 1:n
        r1(ii) = r.^(ii-1);
    end
    
    Ce = toeplitz(r1)*(sigma.^2);
    
    Ce = sparse(Ce);
    
end



end


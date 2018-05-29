function [chainstr] = mkchainstr(ichain)
%[chainstr] = mkchainstr(nchain)
%   Make string of letters unique to the chain number - basically convert
%   ichain to an base 26 number, before converting to a N-character string
%   of letters,

if ichain <= 26
    A = ichain;
    chainstr = char(A+64);
elseif ichain <= 26*27
    B = floor((ichain-1)/26);
    A = ichain - 26*B;
    chainstr = [char(B+64),char(A+64)];
% elseif ichain <= 27*26^2
%     C = floor((ichain-1)/26/27);
%     B = floor(((ichain - 26*26*C)-1)/26);
%     A = ichain - 26*B - 26*26*C;
%     chainstr = [char(C+64),char(B+64),char(A+64)];
end

end


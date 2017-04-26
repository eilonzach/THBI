function [ nchain ] = kchain_addcount( nchain,ptbnorm,par )
%[ nchain ] = kchain_addcount( nchain,ptbnorm,par )
%   function to augment "nchain" counter based on the ptbnorm and the
%   threshold values

if ptbnorm > par.inv.kerneltolmax
    nchain = nchain + 20;
elseif ptbnorm > par.inv.kerneltolmed
    nchain = nchain + 10;
elseif ptbnorm > par.inv.kerneltolmin
    nchain = nchain + 5;  
else
    nchain = nchain + 1;
end


end


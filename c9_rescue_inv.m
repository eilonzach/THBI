function [misfits_perchain,allmodels_perchain] = c9_rescue_inv(par,resdir)
%[misfits_perchain,allmodels_perchain] = c9_rescue_inv(par,resdir)
%   Function to load saved states of chains in the even the inversion is
%   stopped early. 

misfits_perchain = cell(par.inv.nchains,1);
allmodels_perchain = cell(par.inv.nchains,1);

for iii = 1:par.inv.nchains
    chainstr = mkchainstr(iii);
    a = load([resdir,'/',chainstr,'_invState']);

    misfits_perchain{iii} = a.misfits;
    allmodels_perchain{iii} = a.allmodels;
end

end

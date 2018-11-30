function [misfits_perchain,allmodels_perchain,par,trudata,prior] = c9_rescue_inv(resdir)
%[misfits_perchain,allmodels_perchain,par,trudata,prior] = c9_rescue_inv(resdir)
%   Function to load saved states of chains in the even the inversion is
%   stopped early. 

load([resdir,'/par.mat'],'par');

misfits_perchain = cell(par.inv.nchains,1);
allmodels_perchain = cell(par.inv.nchains,1);

for iii = 1:par.inv.nchains
    chainstr = mkchainstr(iii);
    a = load([resdir,'/',chainstr,'_invState']);

    misfits_perchain{iii} = a.misfits;
    allmodels_perchain{iii} = a.allmodels;
end

load([resdir,'/prior.mat'],'prior');
load([resdir,'/trudata_USE.mat'],'trudata');

end


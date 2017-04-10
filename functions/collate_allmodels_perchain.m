function [ allmodels_collated ] = collate_allmodels_perchain( allmodels_perchain,par )
% [ allmodels_collated ] = collate_allmodels_perchain( allmodels_perchain )
% 
%  Function to collate a number of allmodels structures (each occupying a
%  separate cell in allmodels_perchain) into a single structure, with
%  models prior to burn in removed, and with a field that recalls which
%  chain each model came from.

allmodels_collated = [];

for iii = 1:length(allmodels_perchain)
    am = allmodels_perchain{iii};
    if isempty(am), continue; end
    am([am.iter]'<=par.inv.burnin) =  [];
    am = dealto(am,'chain',iii);
    
    allmodels_collated = [allmodels_collated,am];

end

end


function avNmod = count_N_layers_in_mod( allmodels_collated,par )
% avNmod = count_N_layers_in_mod( allmodels_collated )
% 
% sample collated set of models (using only 'bestmods') and layerise a few
% times to calculate the average number of layers in the accepted models

am = allmodels_collated;

bms = find([am.bestmods]');
k = length(bms);
% inds = randperm(length(bms),k);
% bms = bms(inds);

Nmod = zeros(k,1);
for ik = 1:k
    mod = am(bms(ik));
[zlayt,zlayb,Vslay] = ...
    layerise(mod.z,mod.VS,par.forc.mindV,0);
Nmod(ik) = length(zlayt);
end
end


function [ SWwt ] = calc_K_in_model( K,par )
%[ SWwt ] = cald_K_in_model( K )
%   function to calculate the fraction of each kernel that is in the model
%   space - kernels for longer periods will have a significant portion of
%   sensitivity below the model base. We might not want to treat these as
%   equal, as their misfit is heavily contingent on a single parameter -
%   the velocity of our bottom-most mantle node.

Np = length(K);
SWwt = zeros(Np,1);
for ip = 1:Np
    Kf1 = sum(abs(K{ip}.Vsv(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vsv));
    Kf2 = sum(abs(K{ip}.Vsh(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vsh));
    Kf3 = sum(abs(K{ip}.Vpv(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vpv));
    Kf4 = sum(abs(K{ip}.Vph(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vph));
    Kf = [Kf1,Kf2,Kf3,Kf4];
    SWwt(ip) = nansum(Kf)/sum(~isnan(Kf));
end


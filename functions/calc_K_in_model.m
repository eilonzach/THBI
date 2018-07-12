function [ SWwt ] = calc_K_in_model( K,par )
%[ SWwt ] = cald_K_in_model( K )
%   function to calculate the fraction of each kernel that is in the model
%   space - kernels for longer periods will have a significant portion of
%   sensitivity below the model base. We might not want to treat these as
%   equal, as their misfit is heavily contingent on a single parameter -
%   the velocity of our bottom-most mantle node.

Np = length(K);
SWwt = ones(Np,1);

% option for MINEOS-type kernels for phase/group velocities
if isfield(K{1},'Z')
    for ip = 1:Np
        Kf1 = sum(abs(K{ip}.Vsv(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vsv));
        Kf2 = sum(abs(K{ip}.Vsh(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vsh));
        Kf3 = sum(abs(K{ip}.Vpv(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vpv));
        Kf4 = sum(abs(K{ip}.Vph(0.001*K{ip}.Z<par.mod.maxz)))/sum(abs(K{ip}.Vph));
        Kf = [Kf1,Kf2,Kf3,Kf4];
        SWwt(ip) = nansum(Kf)/sum(~isnan(Kf));
    end
% option for Tanimoto-type kernels for HV ratios    
elseif isfield(K{1},'Z1')
    Zav = 0.5*(K{1}.Z1 + K{1}.Z2);
    if max(Zav)>par.mod.maxz % don't even bother with weights if all inside model
        for ip = 1:Np
            Kf1 = sum(abs(K{ip}.Kzh_Vs(0.001*Zav<par.mod.maxz)))/sum(abs(K{ip}.Kzh_Vs));
            Kf2 = sum(abs(K{ip}.Kzh_Vp(0.001*Zav<par.mod.maxz)))/sum(abs(K{ip}.Kzh_Vp));
            Kf3 = sum(abs(K{ip}.Kzh_rho(0.001*Zav<par.mod.maxz)))/sum(abs(K{ip}.Kzh_rho));
            Kf = [Kf1,Kf2,Kf3];
            SWwt(ip) = nansum(Kf)/sum(~isnan(Kf));
        end
    end
end



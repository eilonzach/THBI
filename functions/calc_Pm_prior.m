function Pm_prior = calc_Pm_prior(model,par)
% Pm_prior = calc_Pm_prior(model,par)
% 
% Calculate prior probability of model on the basis of parameters
% 
% NOTE: for now, only certain parameters are checked for prior probability.
% Can expand to more as needed

Pm_prior = 1;

% probabilities are combined in series, multiplicatively (i.e. overall prob
% is product of individual probs - this assumes no covariance between
% parameters.

%% Crustal thickness
try
    Pm_prior = Pm_prior*par.mod.crust.h_pprior(model.zmoh-model.zsed);
catch
    if par.inv.verbose
        fprintf(' >> Model prior prob not being calculated for crust h <<\n');
    end
end

%% Crustal Vp/Vs
try
    Pm_prior = Pm_prior*par.mod.crust.vpvs_pprior(model.vpvs);
catch
    if par.inv.verbose
        fprintf(' >> Model prior prob not being calculated for crust vpvs <<\n');
    end
end



end


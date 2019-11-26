function [ sigma_min,sigma_prior ] = get_sigma_min_prior( dtype,par )
%  [ sigma_min,sigma_prior ] = get_sigma_min_prior( dtype )
% 
% function to get the prior uncertainty information for data dtype

pdt = parse_dtype(dtype);


hypmin = par.mod.data.min_sigma;
hyprior = par.mod.data.prior_sigma;

% min
if strcmp(pdt{1},'BW')
    if strcmp(pdt{4},'lo') % lo takes precedence
        sigma_min   = hypmin.(pdt{1}).(pdt{2}).(pdt{4});
        sigma_prior = hyprior.(pdt{1}).(pdt{2}).(pdt{4});
    elseif strcmp(pdt{3},'cms') % cms next
        sigma_min   = hypmin.(pdt{1}).(pdt{2}).(pdt{3});
        sigma_prior = hyprior.(pdt{1}).(pdt{2}).(pdt{3});
    else % default
        sigma_min   = hypmin.(pdt{1}).(pdt{2}).def;
        sigma_prior = hyprior.(pdt{1}).(pdt{2}).def;
    end
elseif strcmp(pdt{1},'RF')
    if strcmp(pdt{4},'lo') % lo takes precedence
        sigma_min   = hypmin.(pdt{1}).(pdt{2}).(pdt{4});
        sigma_prior = hyprior.(pdt{1}).(pdt{2}).(pdt{4});
    elseif strcmp(pdt{3},{'cms','ccp'}) % cms or ccp next
        sigma_min   = hypmin.(pdt{1}).(pdt{2}).(pdt{3});
        sigma_prior = hyprior.(pdt{1}).(pdt{2}).(pdt{3});
    else % default
        sigma_min   = hypmin.(pdt{1}).(pdt{2}).def;
        sigma_prior = hyprior.(pdt{1}).(pdt{2}).def;
    end
elseif strcmp(pdt{1},'SW')
    % default only, for now (Rayleigh, phV)
    sigma_min   = hypmin.(pdt{1}).(pdt{2}).(pdt{3});
    sigma_prior = hyprior.(pdt{1}).(pdt{2}).(pdt{3});
elseif strcmp(pdt{1},'HKstack')
    % default only, (P)
    sigma_min   = hypmin.(pdt{1}).(pdt{2});
    sigma_prior = hyprior.(pdt{1}).(pdt{2});
end
        

end


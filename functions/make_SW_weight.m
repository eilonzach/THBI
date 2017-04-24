function [ SWwt ] = make_SW_weight( par,Kbase )
%[ SWwt ] = make_SW_weight( par,Kbase)
%   
% function to get surface wave weights if par says to do so, and otherwise
% not weight

if all(par.inv.Kweight == true) % if using default weight by fraction of kernel in model
    SWwt = calc_K_in_model( Kbase.K,par );
    SWwt=SWwt/mean(SWwt);
elseif all(par.inv.Kweight == false) % if no weight
    SWwt = [];
else
    SWwt = par.inv.Kweight; % if not explicitly "false", assume custom wts and use.
end

end


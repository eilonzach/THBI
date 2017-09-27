function [ SWwt ] = make_SW_weight( par,Kbase )
%[ SWwt ] = make_SW_weight( par,Kbase)
%   
% function to get surface wave weights if par says to do so, and otherwise
% not weight

for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id}; pdtyp = parse_dtype(dtype);
    if ~strcmp(pdtyp{1},'SW'),continue; end

    % what should the weight be?
    if all(par.inv.Kweight == true) % if using default weight by fraction of kernel in model
            SWwt.(dtype) = calc_K_in_model( Kbase.(pdtyp{2}).(['K',pdtyp{3}(1:2)]),par );
            SWwt.(dtype) = SWwt.(dtype)/mean(SWwt.(dtype));    
    elseif all(par.inv.Kweight == false) % if no weight
        SWwt.(dtype) = 1;    
    else
        SWwt.(dtype) = par.inv.Kweight; % if not explicitly "false", assume custom wts and use.
    end

end % loop on data types


end


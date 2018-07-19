function [ SWwt ] = make_SW_weight( par,Kbase,trudata )
%[ SWwt ] = make_SW_weight( par,Kbase,trudata)
%   
% function to get surface wave weights if par says to do so, and otherwise
% not weight
SWwt=[];
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id}; pdtyp = parse_dtype(dtype);
    if ~strcmp(pdtyp{1},'SW'),continue; end

    switch pdtyp{2}
    
        case {'Ray','Lov'}
            % what should the weight be?
            if all(par.inv.Kweight == true) % if using default weight by fraction of kernel in model
                    SWwt.(dtype) = calc_K_in_model( Kbase.(pdtyp{2}).(['K',pdtyp{3}(1:2)]),par );
                    SWwt.(dtype) = SWwt.(dtype)/geomean(SWwt.(dtype));    
            elseif all(par.inv.Kweight == false) % if no weight
                SWwt.(dtype) = 1;    
            else
                SWwt.(dtype) = par.inv.Kweight; % if not explicitly "false", assume custom wts and use.
            end
        
        case {'HV'}
            % what should the weight be?
            if all(par.inv.Kweight == true) % if using default weight by fraction of kernel in model
                    SWwt.(dtype) = calc_K_in_model( Kbase.(pdtyp{2}).(['K',pdtyp{3}(1:2)]),par );
                    SWwt.(dtype) = SWwt.(dtype)/geomean(SWwt.(dtype));    
            elseif all(par.inv.Kweight == false) % if no weight
                SWwt.(dtype) = 1;    
            else
                SWwt.(dtype) = par.inv.Kweight; % if not explicitly "false", assume custom wts and use.
            end
            
    end
    % allow for different stds at different periods by scaling weights by
    % indiv values (de-meaned using geometrical mean, so doesn't figure
    % later in the expressions)
    if isempty(trudata.(dtype).sigma),trudata.(dtype).sigma=1;end
    dsig = trudata.(dtype).sigma./geomean(trudata.(dtype).sigma);
    SWwt.(dtype) = SWwt.(dtype).*(dsig.^-2);


end % loop on data types


end


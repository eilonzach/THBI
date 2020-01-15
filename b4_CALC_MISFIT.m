function [ misfit ] = b4_CALC_MISFIT( trudata,predata,par,ifplot,SWwt )
% [ misfit ] = b4_CALC_MISFIT( trudata,predata,par,ifplot,SWwt)
%  
% Calculate cross-convolution misfit between observed data and
% (unit-normalised) predicted data

if nargin < 4 || isempty(ifplot)
    ifplot = 0;
end
if nargin < 5 || isempty(SWwt)
    SWwt = struct([]);
    for id = 1:length(par.inv.datatypes)
        if regexp(par.inv.datatypes{id},'SW')
            SWwt(1).(par.inv.datatypes{id}) = 1;
        end
    end
end

%% Setup misfit 
% this is the SUM OF SQUARES MISFIT VECTOR
misfit = struct('E2',cell2struct(cell(size(par.inv.datatypes)),par.inv.datatypes,2));
    


%% BW misfit
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    pdt = parse_dtype(dtype);
    
    % BODY WAVE 
    % Calculate cross-convolution misfit between observed data and
    % (unit-normalised) predicted data
    if any(strcmp({'BW','RF'},pdt{1}))
        misfit2    = zeros(length(trudata.(dtype)),1);
        stfpow_tru = zeros(length(trudata.(dtype)),1);
        stfpow_pre = zeros(length(trudata.(dtype)),1);
        for itr = 1:length(trudata.(dtype))
            misfit2(itr)    = xconv_misfit(trudata.(dtype)(itr).PSV(:,1),...
                                               trudata.(dtype)(itr).PSV(:,2),...
                                               predata.(dtype)(itr).PSV(:,1),...
                                               predata.(dtype)(itr).PSV(:,2));
            stfpow_tru(itr) = norm(trudata.(dtype)(itr).PSV(:,1)) +...
                                 norm(trudata.(dtype)(itr).PSV(:,2));
            stfpow_pre(itr) = norm(predata.(dtype)(itr).PSV(:,1)) +...
                                 norm(predata.(dtype)(itr).PSV(:,2));
        end
        misfit.E2.(dtype) = misfit2./stfpow_tru./stfpow_pre;      % SUM OF SQUARED MISFITS, NORMALISED

    % SURFACE WAVE
    elseif strcmp(pdt{1},'SW')
        e = (trudata.(dtype).(pdt{3}) - predata.(dtype).(pdt{3}));
        misfit.E2.(dtype) = e'*(SWwt.(dtype).*e); % SUM OF SQUARED MISFITS, NORMALISED

    % HKstack
    elseif strcmp(pdt{1},'HKstack')
        % a is deliberate dummy - not sure how to calc error here
        misfit.E2.(dtype) = 1./predata.(dtype).E_by_Emax; % pick out error value

    end % end on data type
end

if ifplot
    plot_TRUvsPRE( trudata,predata)
end

end

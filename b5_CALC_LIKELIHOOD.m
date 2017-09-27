function [ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,data,datahparm,par )
%  [ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,data,datahparm,par )
% 
%   Function to calculate the likelihood function (the inverse exponential
%   of the error-weighted misfit function)
% 
% INPUTS:
%  misfit = misfit structure with squared misfit of each data type
%  data   = true data structure, including estimated sigma of each data
%           type
%  par    = input parameters, including what data types to solve for

% FOR each datatype, do the following:
%  1) get data error
%  2) Number of datapoints
%  3) chi2 of each data type
%  4) save chi2 into misfit structure
%  5) save rms into misfit structure
%  6) calc log_likelihood
% 	   ==> likelihood of each data type
%         sq2p = sqrt(2*pi);
%             we will ignore a constant factor of 
%             -0.5*N_i*log(sq2p)        for each of the data types, "i"
%             because it never changes, and so does not affect results


%% Calc log Likelihood
M = zeros(length(par.inv.datatypes),1);
logL = zeros(length(par.inv.datatypes),1);
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    pdt = parse_dtype(dtype);
    
    sig = datahparm.(['sig_',dtype]);
    M(id) = length(data.(dtype));
    % calculate N - # of periods if SW, # of degrees of freedom if BW
	N = zeros(M(id),1);
    for itr = 1:length(data.(dtype))
        if strcmp(pdt{1},'SW') 
            N(itr) = length(data.(dtype)(itr).periods);
        elseif strcmp(pdt{1},'BW') 
            N(itr) = mean([scdofcalc(data.(dtype)(itr).PSV(:,1)),scdofcalc(data.(dtype)(itr).PSV(:,2))]);
        end
    end
    
    chi2 = misfit.E2.(dtype)./sig./sig;
    rms  = sqrt(misfit.E2.(dtype)./N);
    logL(id) = sum(-N.*log(sig) - 0.5*chi2)./M(id);
    
    misfit.chi2.(dtype) = chi2;
    misfit.rms.(dtype) = rms;
end

%% collate

% summed misfit
misfit.chi2sum = 0;
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    misfit.chi2sum = misfit.chi2sum + sum(misfit.chi2.(dtype))./M(id);
end

% summed log likelihood

log_likelihood = sum(logL);



end


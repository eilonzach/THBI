function   [misfits,allmodels,savedat] = b9_SAVE_RESULT(iter,log_likelihood,misfit,model,Pm_prior,misfits,allmodels,predata,savedat,time0)
%  [misfits,allmodels,savedat] = b9_SAVE_RESULT(iter,log_likelihood,misfit,model,misfits,allmodels,predata,time0)
% 
% Function to store model into results structure and to store the misfit
% for future reference

%% MISFITS
% N = misfits.Nstored+1; 

% choose the next N as the save number that corresponds to one more than
% the numnber saved for iterations less than the current iteration. That
% way if we reset the model of the MCMC, it will back-track and write over
% old saves, too. 
N = sum(misfits.iter < iter & misfits.iter~=0) + 1;

likelihood = exp(log_likelihood);
misfits.lastlogL = log_likelihood;
misfits.lastL = likelihood;
misfits.globmaxL = max(misfits.globmaxL,likelihood);

misfits.logLike(N,1) = log_likelihood;
misfits.Like(N,1) = likelihood;

misfits.chi2sum(N,1)    = misfit.chi2sum;

misfits.chi2(N,1) = misfit.chi2; % chi2 is the chi-squared misfit for each data type, accounting for data error, i.e. E2/sig/sig
misfits.rms(N,1) = misfit.rms; % rms is the root mean squared error for each dat type, i.e. sqrt(E2/N)
misfits.E2(N,1) = misfit.E2; % E2 is the normalised, weighted, sum of squared errors

misfits.time(N,1) = (now-time0)*86400; % in seconds
misfits.iter(N,1) = iter;
misfits.Nstored = N;

%% MODEL
% N = allmodels(1).Nstored+1;
fns = fieldnames(model);
for ifn = 1:length(fns)
    allmodels(N).(fns{ifn}) = model.(fns{ifn});
end
allmodels(N).iter = iter;
allmodels(N).Pm_prior = Pm_prior;
allmodels(1).Nstored = N;

%% DATA
dtypes = fieldnames(predata);
for id = 1:length(dtypes)
    dtype = dtypes{id};
    pdt = parse_dtype(dtype);
    if strcmp(pdt{1},'SW')
        inds = 1:length(predata.(dtype).periods);         
        savedat.(dtype).(pdt{3})(N,inds) = predata.(dtype).(pdt{3});
        savedat.(dtype).periods(N,inds) = predata.(dtype).periods;
    elseif strcmp(pdt{1},'BW')
        for jj = 1:length(predata.(dtype))
            inds = 1:predata.(dtype)(jj).nsamp;
            savedat.(dtype)(jj,1).P(N,inds) = predata.(dtype)(jj).PSV(:,1);
            savedat.(dtype)(jj,1).SV(N,inds) = predata.(dtype)(jj).PSV(:,2);
            savedat.(dtype)(jj,1).tt(N,inds) = predata.(dtype)(jj).tt;
        end
    end
end

end


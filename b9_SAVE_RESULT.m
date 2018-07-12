function   [misfits,allmodels,savedat] = b9_SAVE_RESULT(iter,log_likelihood,misfit,model,misfits,allmodels,predata,savedat)
%  [misfits,allmodels,savedat] = b9_SAVE_RESULT(iter,log_likelihood,misfit,model,misfits,allmodels,predata)
% 
% Function to store model into results structure and to store the misfit
% for future reference

%% MISFITS
N = misfits.Nstored+1;
likelihood = exp(log_likelihood);
misfits.lastlogL = log_likelihood;
misfits.lastL = likelihood;
misfits.globmaxL = max(misfits.globmaxL,likelihood);

misfits.iter(N,1) = iter;
misfits.logLike(N,1) = log_likelihood;
misfits.Like(N,1) = likelihood;

misfits.chi2sum(N,1)    = misfit.chi2sum;

misfits.chi2(N,1) = misfit.chi2;
misfits.rms(N,1) = misfit.rms;
misfits.E2(N,1) = misfit.E2;

misfits.Nstored = N;

%% MODEL
N = allmodels(1).Nstored+1;
fns = fieldnames(model);
for ifn = 1:length(fns)
    allmodels(N).(fns{ifn}) = model.(fns{ifn});
end
allmodels(N).iter = iter;
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


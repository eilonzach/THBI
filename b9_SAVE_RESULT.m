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

misfits.chi2(N,1)    = misfit.chi2;

datatypes = fieldnames(predata);
if any(strcmp(datatypes,'SW')) && isfield(misfit,'chi2_SW')
misfits.chi2_SW(N,1) = misfit.chi2_SW;
misfits.norm_SW(N,1) = misfit.SW;
misfits.rms_SW(N,1) = misfit.rms_SW;
end
if any(strcmp(datatypes,'PsRF')) && isfield(misfit,'chi2_ps')
misfits.chi2_ps(N,:) = misfit.chi2_ps;
misfits.norm_ps(N,:) = misfit.PsRF;
misfits.rms_ps(N,:) = misfit.rms_ps;
end
if any(strcmp(datatypes,'SpRF')) && isfield(misfit,'chi2_sp')
misfits.chi2_sp(N,:) = misfit.chi2_sp;
misfits.norm_sp(N,:) = misfit.SpRF;
misfits.rms_sp(N,:) = misfit.rms_sp;
end
if any(strcmp(datatypes,'PsRF_lo')) && isfield(misfit,'chi2_ps_lo')
misfits.chi2_ps_lo(N,:) = misfit.chi2_ps_lo;
misfits.norm_ps_lo(N,:) = misfit.PsRF_lo;
misfits.rms_ps_lo(N,:) = misfit.rms_ps_lo;
end
if any(strcmp(datatypes,'SpRF_lo')) && isfield(misfit,'chi2_sp_lo')
misfits.chi2_sp_lo(N,:) = misfit.chi2_sp_lo;
misfits.norm_sp_lo(N,:) = misfit.SpRF_lo;
misfits.rms_sp_lo(N,:) = misfit.rms_sp_lo;
end

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
    if regexp(dtype,'SW')
        inds = 1:length(predata.(dtype).periods);         
        savedat.(dtype).phV(N,inds) = predata.SW.phV;
        savedat.(dtype).periods(N,inds) = predata.SW.periods;
    elseif regexp(dtype,'RF')
        if regexp(dtype,'_lo'), continue; end
        for jj = 1:length(predata.(dtype))
            inds = 1:predata.(dtype)(jj).nsamp;
            savedat.(dtype)(jj,1).Z(N,inds) = predata.(dtype)(jj).ZRT(:,1);
            savedat.(dtype)(jj,1).R(N,inds) = predata.(dtype)(jj).ZRT(:,2);
            savedat.(dtype)(jj,1).T(N,inds) = predata.(dtype)(jj).ZRT(:,3);
            savedat.(dtype)(jj,1).tt(N,inds) = predata.(dtype)(jj).tt;
        end
    end
end

end


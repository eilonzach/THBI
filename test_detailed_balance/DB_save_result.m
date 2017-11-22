function   [allmodels] = DB_save_result(iter,log_likelihood,model,allmodels)
%  [allmodels] = DB_save_result(iter,log_likelihood,model,allmodels)
% 
% Function to store model into results structure for future reference

%% MODEL
N = allmodels(1).Nstored+1;
fns = fieldnames(model);
for ifn = 1:length(fns)
    allmodels(N).(fns{ifn}) = model.(fns{ifn});
end
allmodels(N).iter = iter;
allmodels(N).logL = log_likelihood;
allmodels(1).Nstored = N;

end


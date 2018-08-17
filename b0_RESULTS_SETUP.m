function [ misfits,allmodels,savedat,log_likelihood ] = b0_RESULTS_SETUP(par )
% [ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par )
% 
%  Function to set up the results structures for the inversion

%% OUTPUT
o=[];
Nsave = ceil(par.inv.niter/par.inv.saveperN);
nn = nan(Nsave,1);
oo = zeros(Nsave,1);

% misfits 
misfits = struct('globmaxL',0,'lastL',0,'lastlogL',-Inf,... % minimum global errror, most recent error
                 'logLike',nn,'Like',nn,... % likelihoods at each stage in the inversion
                 'chi2sum',nn,...  % chi2sum is the chi-squared misfit summed across all data types
                 'time',nn,... % time since inversion started (in seconds)
                 'iter',0,... % iteration number
                 'Nstored',0);
% Will also store (but done inside b9_SAVE_RESULT)
%       chi2 = the chi-squared misfit for each data type, accounting for data error, i.e. E2/sig/sig
%       rms  = the root mean squared error for each dat type, i.e. sqrt(E2/N)
%       E2   = the normalised, weighted, sum of squared errors
             
             
% log likelihood
log_likelihood = -Inf;

% allmodels
allmodels = struct('z',o,'z0',o,'iter',o,...
               'VS',o,'VP',o,'rho',o,...
               'Nz',o,'zsed',o,'zmoh',o,...
               'fdVSsed',o,'fdVSmoh',o,...
               'Sanis',o,'Panis',o,...
               'sedmparm',struct('h',o,'vstop',o,'vsbot',o),...
               'crustmparm',struct('h',o,'Nkn',o,'VS_kn',o,'splines',o),...
               'mantmparm',struct('Nkn',o,'VS_kn',o,'splines',o),...
               'M',nan,'Nstored',0);
allmodels(Nsave)=allmodels;

%save data
dtypes = par.inv.datatypes;
savedat=cell2struct(cell(1,length(dtypes)),dtypes,2);
for id = 1:length(dtypes)
    dtype = dtypes{id}; pdtyp = parse_dtype(dtype); 
    if regexp(dtypes{id},'SW')
        savedat.(dtypes{id})=struct(pdtyp{3},oo,'periods',oo);
    elseif regexp(dtypes{id},'RF')
%         if regexp(dtypes{id},'_lo'), continue; end
        savedat.(dtypes{id})=struct('P',oo,'SV',oo,'tt',oo);
    end
end

end


function [ misfits,allmodels,savedat,log_likelihood ] = b0_RESULTS_SETUP(par )
% [ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par )
% 
%  Function to set up the results structures for the inversion

%% OUTPUT
o=[];
nn = nan(par.inv.niter,1);
oo = zeros(par.inv.niter,1);

% misfits 
misfits = struct('globmaxL',0,'lastL',0,'lastlogL',-Inf,... % minimum global errror, most recent error
                 'iter',0,... % iteration number
                 'chi2sum',nn,...  % chi2 is the chi-squared misfit, accounting for data error
                 'Nstored',0);
                 % norm is the sum of the squared error
                 
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
allmodels(par.inv.niter)=allmodels;

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


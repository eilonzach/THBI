function posterior = c2_BUILD_POSTERIOR(allmodels,par,zatdep)
% profile clear
% profile on

o = [];

posterior = struct('Niter',par.inv.niter*par.inv.nchains,'Nstored',0,...
               'zsed',o,'zmoh',o,...
               'kcrust',o,'kmantle',o,...
               'VSsedtop',o,'VSsedbot',o,'VScrusttop',o,'VScrustbot',o,'VSmanttop',o,...
               'VSmantle',repmat(o,1,6),'zmantle',o,'NVGzwa',o,...
               'fdVSsed',o,'fdVSmoh',o,'vpvs',o,'xicrust',o,'ximant',o,...
               'datahparm',allmodels(1).datahparm);
if nargin<3 || isempty(zatdep)
    posterior.zatdep = linspace(par.mod.sed.hmax+par.mod.crust.hmin+0.1,par.mod.maxz,50)';
else
    posterior.zatdep = zatdep;
end

dtypes = fieldnames(allmodels(1).datahparm);

if iscell(allmodels)
    nchains = length(allmodels);
else
    nchains = 1;
end


for iii = 1:nchains
    if nchains>1, am = allmodels{iii}; else am = allmodels; end
    am([am.bestmods]'==false) =  [];

    % get indices of posterior row to put result in - accounting for
    % previous entries from other chains
    ind0 = posterior.Nstored+1;
    ind1 = posterior.Nstored+length(am);
    inds = ind0:ind1;
    
    % put into 1-D prior struct
    posterior.Nstored = ind1;
    posterior.zsed(inds,1) = [am.zsed]';
    posterior.zmoh(inds,1) = [am.zmoh]';
    posterior.fdVSsed(inds,1) = [am.fdVSsed]';
    posterior.fdVSmoh(inds,1) = [am.fdVSmoh]';
    posterior.zatdep = round(posterior.zatdep);
    for ii = 1:length(am)
        posterior.kcrust(inds(ii),1) = am(ii).crustmparm.Nkn;
        posterior.kmantle(inds(ii),1) = am(ii).mantmparm.Nkn;
        posterior.VSsedtop(inds(ii),1) = am(ii).sedmparm.VS(1);
        posterior.VSsedbot(inds(ii),1) = am(ii).sedmparm.VS(2);
        posterior.VScrusttop(inds(ii),1) = am(ii).crustmparm.VS_sp(1);
        posterior.VScrustbot(inds(ii),1) = am(ii).crustmparm.VS_sp(end);
        posterior.VSmanttop(inds(ii),1) = am(ii).mantmparm.VS_sp(1);
        posterior.VSmantle(inds(ii),:) = linterp(am(ii).z,am(ii).VS,posterior.zatdep);
        posterior.vpvs(inds(ii),1) = am(ii).crustmparm.vpvs;
        posterior.xicrust(inds(ii),1) = am(ii).crustmparm.xi;
        posterior.ximant(inds(ii),1) = am(ii).mantmparm.xi;
        for id = 1:length(dtypes)
            posterior.datahparm.(dtypes{id})(:,inds(ii)) = am(ii).datahparm.(dtypes{id})';
        end
        posterior.ximant(inds(ii),1) = am(ii).mantmparm.xi;
        posterior.ximant(inds(ii),1) = am(ii).mantmparm.xi;
        
        % nvg
        [zwa(1),zwa(2),zwa(3)] =  model_NVG_info(am(ii));
        posterior.NVGzwa(inds(ii),:) = zwa;
        
            
    end    
    
end

% flip the datahparm sigma values to be columns
for id = 1:length(dtypes)
    posterior.datahparm.(dtypes{id}) = posterior.datahparm.(dtypes{id})';
end


end
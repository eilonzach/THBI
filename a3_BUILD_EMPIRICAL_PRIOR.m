function empprior = a3_BUILD_EMPIRICAL_PRIOR(par,Niter,Nchains,zatdep)
% empprior = a3_BUILD_EMPIRICAL_PRIOR(par,Niter,Nchains,zatdep)

if nargin <2 || isempty(Niter) 
    Niter = 1e4;
end
if nargin <3 || isempty(Nchains) 
    Nchains = 10;
end
if nargin<4 || isempty(zatdep)
    zatdep = linspace(par.mod.sed.hmax+par.mod.crust.hmin+0.1,par.mod.maxz,50)';
end
           
% zatdep = [5:5:par.mod.maxz]';
Nsave = floor(par.inv.niter/par.inv.saveperN);
o = nan(Nsave,1);


   
%% START DIFFERENT MARKOV CHAINS IN PARALLEL
%% ========================================================================
%% ========================================================================
fprintf('\n =========== STARTING PARALLEL CHAINS ===========\n')
%% ========================================================================
%% ========================================================================
t = now;
parfor iii = 1:Nchains 

chainstr = mkchainstr(iii);

%% Initiate model
ifpass = 0;
% only let starting model out of the loop if it passes conditions
while ifpass==0 % first make sure the starting model satisfies conditions
    model0 = b1_INITIATE_MODEL(par);
    model = model0;
    ifpass = a1_TEST_CONDITIONS( model, par );
end % now we have a starting model!


dtypes = fieldnames(model(1).datahparm);
empprior = struct('Niter',Niter*Nchains,'Nstored',0,...
               'zsed',o,'zmoh',o,...
               'kcrust',o,'kmantle',o,...
               'VSsedtop',o,'VSsedbot',o,'VScrusttop',o,'VScrustbot',o,'VSmanttop',o,...
               'VSmantle',nan(Nsave,length(zatdep)),'zatdep',zatdep,...
               'Zkn_crust',nan(Nsave,par.mod.crust.kmax),'Zkn_mantle',nan(Nsave,par.mod.mantle.kmax),...
               'fdVSmoh',o,'vpvs',o,'cxi',o,'mxi',o,'NVGzwa',o,...
               'datahparm',nan(Nsave,length(dtypes)));

%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
ptb = cell({});
acc = zeros(par.inv.niter,1);
ifaccept=false; 
% reset_likelihood;
lastlogL = -Inf; 
% not parfor
fprintf('\n =========== STARTING ITERATIONS %s ===========\n',chainstr)
for ii = 1:Niter
    if rem(ii,1000)==0, fprintf('Iteration %s%.0f\n',chainstr,ii); end
    temp = 1;
    
%% ===========================  PERTURB  ===========================  
    if ii==1
        model1 = model; % don't perturb on first run
        p_bd = 1;
        ptb{ii,1} = 'start';
    else
        [model1,ptb{ii,1},p_bd] = b2_PERTURB_MODEL(model,par,temp);
        ifpass = a1_TEST_CONDITIONS( model1, par, par.inv.verbose );
        if ifpass==0 && model1.mantmparm.knots(end-1)<par.mod.maxkz
        end
    end
 
%% =====================    ===  ACCEPTANCE CRITERION  ========================
    [ ifaccept ] = b6_IFACCEPT( 0,lastlogL,temp,p_bd*ifpass );

        
%% ========================  IF ACCEPT ==> STORE  =========================
    if ifaccept 
        model = model1;
        lastlogL=0;
        acc(ii)=1;
    else
        if par.inv.verbose, fprintf('    > fail\n'); end
        acc(ii) = 0;
    end
    
    
    % SAVE into 1-D empprior struct
    if mod(ii,par.inv.saveperN)==0
    istor = empprior.Nstored+1;
    empprior.Nstored=istor;
    empprior.ifaccept(istor,1) = ifaccept;
    empprior.zsed(istor,1) = model.zsed;
    empprior.zmoh(istor,1) = model.zmoh;
    empprior.fdVSmoh(istor,1) = model.fdVSmoh;
    empprior.kcrust(istor,1) = model.crustmparm.Nkn;
    empprior.kmantle(istor,1) = model.mantmparm.Nkn;
    empprior.VSsedtop(istor,1) = model.sedmparm.VS(1);
    empprior.VSsedbot(istor,1) = model.sedmparm.VS(2);
    empprior.VScrusttop(istor,1) = model.crustmparm.VS_sp(1);
    empprior.VScrustbot(istor,1) = model.crustmparm.VS_sp(end);
    empprior.VSmanttop(istor,1) = model.mantmparm.VS_sp(1);
    empprior.VSmantle(istor,:) = linterp(model.z,model.VS,empprior.zatdep);
    empprior.Zkn_crust(istor,1:model.crustmparm.Nkn-2) = model.crustmparm.knots(2:end-1);
    empprior.Zkn_mantle(istor,1:model.mantmparm.Nkn-2) = model.mantmparm.knots(2:end-1);
    empprior.vpvs(istor,1) = model.crustmparm.vpvs;
    empprior.cxi(istor,1) = model.crustmparm.xi;
    empprior.mxi(istor,1) = model.mantmparm.xi;
    for id = 1:length(dtypes)
        empprior.datahparm(istor,id) = model.datahparm.(dtypes{id});
    end
    
    % nvg
    [empprior.NVGzwa(istor,1),...
        empprior.NVGzwa(istor,2),...
            empprior.NVGzwa(istor,3)] =  model_NVG_info(model);
        
    
    
    end

end % on iterations
%% -------------------------- End iteration  ------------------------------
fprintf('\n =========== ENDING ITERATIONS %s ===========\n',char(64+iii))
emppriors{iii} = empprior;
ptbs{iii} = ptb;

end % parfor loop
%% ========================================================================
%% ========================================================================
%% ----------------------- End loop on chains  ----------------------------
%% ========================================================================
%% ========================================================================
fprintf('Duration of entire run: %.0f s\n',(now-t)*86400)
delete(gcp('nocreate'))

empprior = emppriors{1};
ptb = ptbs{1};
fns = fieldnames(empprior);
for iii = 2:Nchains 
    for jj = 1:length(fns)
        if isscalar(empprior.(fns{jj}))
            empprior.(fns{jj}) = empprior.(fns{jj}) + emppriors{iii}.(fns{jj});
        elseif strcmp(fns{jj},'zatdep')
            continue;
        else
            empprior.(fns{jj}) = [empprior.(fns{jj}) ; emppriors{iii}.(fns{jj})];
        end
    end
    ptb = {ptb;ptbs{iii}};
end

end
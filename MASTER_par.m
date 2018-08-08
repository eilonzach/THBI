clear all
close all


projname = 'SYNTHETICS'; % SYNTHETICS, LAB_tests, or WYOMING, for now
sta = 'EYMN';
nwk = 'US';
gc = [69,59,40,38,36,66]; % will search for gcarcs +/-3 of this value;
% baz = 315;

notes = [...
    'Synthetic test using spline-based model. '...
    'Using all data types. '...
    'Changed the way the spline kernels are re-calculated after N iterations, so as to prevent big jumps in likelihood and therefore hanging chains. '...
    'Several previous runs have shut down computer for reasons I do not understand. '...
    'Execute for few iterations to try to see what breaks, but in a short run time... '...
]; 

%% ------------------------- START ------------------------- 
global projdir THBIpath TRUEmodel
THBIpath = '/Users/zeilon/Documents/MATLAB/BayesianJointInv';
projdir = [THBIpath,'/',projname,'/'];
cd(projdir);

run([THBIpath,'/a0_STARTUP_BAYES']);
run('project_details');


%% PARMS
run parms/bayes_inv_parms
if strcmp(projname,'SYNTHETICS')
    if isfield(par.synth,'noisetype') && strcmp(par.synth.noisetype,'real'), sta=['SYNTH_',sta]; else sta = 'SYNTH'; end
    par.sta = sta; 
    par.nwk = '--';

	% noise details, if "real"
	noisesta = 'RSSD';
	noisenwk = 'IU';
	noisegcarcs = [73,38];
	noiseshape = 'real'; % 'white' or 'real'
	noiseup = 0.5; % factor to increase real noise

    ifsavedat = false;

elseif strcmp(projname,'LAB_tests')
	zsed = 0;
	zmoh = 45;
	zlab = 130;
	wlab = 10;
	flab = 0.05;
    par.synth.model = struct('zsed',zsed,'zmoh',zmoh,'zlab',zlab,'wlab',wlab,'flab',flab);
	dtps = {'BW_Ps','BW_Sp','BW_Sp_lo','BW_Ps_lo','SW_Ray_phV','SW_Lov_phV','SW_HV'};     

	% noise details, if "real"
	noisesta = 'RSSD';
	noisenwk = 'IU';
	noisegcarcs = [73,38];
	noiseshape = 'real'; % 'white' or 'real'
	noiseup = 0.5; % factor to increase real noise

	% naming convention
	dtpstr='_';
	if any(strcmp(dtps,'BW_Ps')), dtpstr=[dtpstr,'Ps']; end
	if any(strcmp(dtps,'BW_Ps_lo')), dtpstr=[dtpstr,'Pslo']; end
	if any(strcmp(dtps,'BW_Ps_cms')), dtpstr=[dtpstr,'Pscms']; end
	if any(strcmp(dtps,'BW_Sp')), dtpstr=[dtpstr,'Sp']; end
	if any(strcmp(dtps,'BW_Sp_lo')), dtpstr=[dtpstr,'Splo']; end
	if any(strcmp(dtps,'SW_Ray_phV')), dtpstr=[dtpstr,'Ray']; end
	if any(strcmp(dtps,'SW_Lov_phV')), dtpstr=[dtpstr,'Lov']; end
	if any(strcmp(dtps,'SW_HV')), dtpstr=[dtpstr,'HV']; end

	ifsavedat = true;

	sta = ['LAB_s',num2str(zsed),'_m',num2str(zmoh),'_z',num2str(zlab),'_w',num2str(wlab),'_f',num2str(100*flab),dtpstr];
    par.sta = sta; par.nwk = 'LAB_test';
else
    par.sta = sta;
    par.nwk = nwk;
    par.gc = gc;
    datN = 20;
end

par.synth.noise_sta_deets = struct('datadir',['/Volumes/data/THBI/US/STAsinv/',noisesta,'_dat20/'],...
                         'sta',noisesta,'nwk',noisenwk,'gc',noisegcarcs,'noiseup',noiseup,'noiseshape',noiseshape);


par.inv.verbose=false;

%% get saving things ready
STAMP=[sta,datestr(now,'_yyyymmddHHMM_pll')]; par.STAMP = STAMP;
avardir = sprintf('%s/%s_dat%.0f/',proj.STAinversions,sta,datN);
resdir = [avardir,STAMP];
mkdir(resdir);
fid = fopen([resdir,'/notes.txt'],'w'); fprintf(fid,notes); fclose(fid);

par_ORIG = par;
% par came from above script
save([resdir,'/par'],'par');
copyfile('parms/bayes_inv_parms.m',resdir);

for id = 1:length(par.inv.datatypes)
    allpdytp(id,:)=parse_dtype(par.inv.datatypes{id});
end

%% ========================  LOAD + PREP DATA  ========================  
trudata = a2_LOAD_DATA(par,projname,resdir);


%% ===========================  PRIOR  ===========================  
zatdep = [5:5:par.mod.maxz]';
% fprintf('  > Building prior distribution from %.0f runs\n',max([par.inv.niter,1e5]))
% prior = a3_BUILD_EMPIRICAL_PRIOR(par,max([par.inv.niter,1e5]),14,zatdep);
% plot_MODEL_SUMMARY(prior,1,[resdir,'/prior_fig.pdf']);
% save([resdir,'/prior'],'prior','par');

%% ---------------------------- INITIATE ----------------------------
%% ---------------------------- INITIATE ----------------------------
%% ---------------------------- INITIATE ----------------------------
profile clear
profile on
delete(gcp('nocreate'));

TD = parallel.pool.Constant(trudata);
% TD(1).Value = trudata;

%% START DIFFERENT MARKOV CHAINS IN PARALLEL
% nchains = ;
% ifpass = zeros(nchains,1);
model0_perchain = cell(par.inv.nchains,1);
misfits_perchain = cell(par.inv.nchains,1);
allmodels_perchain = cell(par.inv.nchains,1);
SWs_perchain = cell(par.inv.nchains,1);

%% ========================================================================
%% ========================================================================
fprintf('\n =========== STARTING PARALLEL CHAINS ===========\n')
%% ========================================================================
%% ========================================================================
t = now;
% mkdir([resdir,'/chainout']);
parfor iii = 1:par.inv.nchains 
chainstr = mkchainstr(iii);


%% Fail-safe to restart chain if there's a succession of failures
fail_chain=20;
while fail_chain>=20

%% Prep posterior structure
[ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par); %#ok<PFTUS>

%% Initiate model
ifpass = 0;
% only let starting model out of the loop if it passes conditions
while ifpass==0
    while ifpass==0 % first make sure the starting model satisfies conditions
        model0 = b1_INITIATE_MODEL(par);
        model = model0;
        model = expand_dathparms_to_data( model,TD.Value,par );
        ifpass = a1_TEST_CONDITIONS( model, par );
    end

    %% starting model kernel
    fprintf('\nCreating starting kernels %s\n',chainstr)
    try
        [Kbase] = make_allkernels(model,[],TD.Value,['start',chainstr],par);
    catch
        ifpass = false; continue; 
    end
    
    model0_perchain{iii} = model0;
    
end % now we have a starting model!

%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
ptb = cell({});
nchain = 0;
fail_chain = 0;
ifaccept=true;
if isfield(TD.Value,'SW_Ray')
    preSW = zeros(length(TD.Value.SW_Ray.periods),ceil(par.inv.niter./par.inv.saveperN));
end
% reset_likelihood;
log_likelihood = -Inf;
predata=[]; predat_save = []; misfit = [];
% not parfor
fprintf('\n =========== STARTING ITERATIONS %s ===========\n',chainstr)
for ii = 1:par.inv.niter
    
%% SAVE model every saveperN
if mod(ii,par.inv.saveperN)==0 && log_likelihood ~= -Inf
    [misfits,allmodels,savedat] = b9_SAVE_RESULT(ii,log_likelihood,misfit,model,misfits,allmodels,predat_save,savedat);
    if isfield(TD.Value,'SW_Ray')
        preSW(:,misfits.Nstored) = predata.SW_Ray.phV;
    end
end

%% SAVE inv state every 2500 iterations
if rem(ii,2500)==0
    save_inv_state(resdir,chainstr,allmodels,misfits)
end
    
try
    if rem(ii,4*par.inv.saveperN)==0 || ii==1, fprintf('Iteration %s%.0f\n',chainstr,ii); end
    if par.inv.verbose, pause(0.05); end
    ifaccept=false;
    ifpass = false;
    newK = false;
    if fail_chain>19
        % if not enough saved in this chain, abort and restart
        if (ii - par.inv.burnin)/par.inv.saveperN < 200
            break
        % if enough saved in chain, abort and keep the incomplete chain
        else 
            fail_chain = -fail_chain; break
        end
    end
    
    % temperature - for perturbation scaling and likelihood increase
    temp = (par.inv.tempmax-1)*erfc(2*(ii-1)./par.inv.cooloff) + 1;
%     if round_level(temp,0.01)>1
%         if par.inv.verbose, fprintf('TEMPERATURE = %.2f\n',temp); end
%     end
    
    while ifpass == false % only keep calculating if model passes (otherwise save and move on)

%% ===========================  PERTURB  ===========================  
    if ii==1
        model1 = model; % don't perturb on first run
        ptbnorm = 0;
        ifpass = 1;
        p_bd = 1;
        log_likelihood1 = -Inf;
        ptb{ii,1} = 'start';
    else
		[model1,ptb{ii,1},p_bd] = b2_PERTURB_MODEL(model,par,temp);
		ifpass = a1_TEST_CONDITIONS( model1, par, par.inv.verbose  );
		if ~ifpass, if par.inv.verbose, fprintf('  nope\n'); end; break; end
		
		[ modptb ] = calc_Vperturbation( Kbase.modelk,model1);
		ptbnorm = 0.5*(norm(modptb.dvsv) + norm(modptb.dvsh)) + norm(modptb.dvpv);
    end

    if isempty(gcp('nocreate'))
%         if exist('h','var'), delete(h); end
        if par.inv.verbose 
        figure(85);clf,set(gcf,'pos',[1294 4 622 529])
        subplot(1,2,1), hold on, 
        plot(TRUEmodel.VS,TRUEmodel.z,'k','linewidth',1); set(gca,'ydir','reverse')
        plot(model.VS,model.z,'r','linewidth',1.5); set(gca,'ydir','reverse')
        plot(model1.VS,model1.z,'b--','linewidth',1.5); set(gca,'ydir','reverse')
        subplot(1,2,2), hold on, 
        plot(TRUEmodel.VP,TRUEmodel.z,'k','linewidth',1); set(gca,'ydir','reverse')
        plot(model.VP,model.z,'r','linewidth',1.5); set(gca,'ydir','reverse')
        plot(model1.VP,model1.z,'b--','linewidth',1.5); set(gca,'ydir','reverse')
        pause(0.01)
        end
    end    
    nchain  = kchain_addcount( nchain,ptbnorm,par );
    
%% ===========================  FORWARD MODEL  ===========================
	% don't re-calc if the only thing perturbed is the error
    if ~strcmp('sigma',ptb{ii}(end-4:end)) || isempty(predata)
        % make random run ID (to avoid overwrites in parfor)
		ID = [chainstr,num2str(round(1e9*(now-t))),num2str(randi(9)),num2str(randi(9))];

        try
            predata = b3_FORWARD_MODEL_BW( model1,par,TD.Value,ID,0 );
            predata = b3_FORWARD_MODEL_SW_kernel( model1,Kbase,par,predata );
        catch
            fail_chain=fail_chain+1;
            fprintf('Forward model error, failchain %.0f\n',fail_chain);  break;
        end
        
        % continue if any Sp or PS inhomogeneous or nan or weird output
        if ifforwardfail(predata,par)
            fail_chain=fail_chain+1; 
            fprintf('Forward model error, failchain %.0f\n',fail_chain);  break;
        end
        
        predata0=predata;
        
        for idt = 1:length(par.inv.datatypes)
            predata = predat_process( predata,par.inv.datatypes{idt},par);
        end
        
		% Explicitly use mineos + Tanimoto scripts if ptb is too large
		if par.inv.verbose, fprintf('    Perturbation %.2f\n',ptbnorm); end
        if ptbnorm/par.inv.kerneltolmax > random('unif',par.inv.kerneltolmed/par.inv.kerneltolmax,1,1) % control chance of going to MINEOS
            newK = true;
        end
%         % ALSO use MINEOS if chain too long
        if nchain > 2*par.inv.maxnkchain
            fprintf('Chain %s, iter %.0f, nchain=%.0f\n',chainstr,ii,nchain);
        end
            
            
        if newK
            [ predata,HVK_new ] = b3_FORWARD_MODEL_SW_precise( model1,par,predata,ID );
		else 
			Ktry = [];
        end
    end % only redo data if model has changed 

%      plot_TRUvsPRE_old(TD.Value,predata)]

    % continue if any Sp or PS inhomogeneous or nan or weird output
    if ifforwardfail(predata,par)
        fail_chain=fail_chain+1; ifpass=0; 
        fprintf('Forward model error, failchain %.0f\n',fail_chain);  break;
    else
        fail_chain = 0; 
    end

%% =========================  CALCULATE MISFIT  ===========================
    
    % SW weights, if applicable 
    [ SWwt ] = make_SW_weight( par,Kbase,TD.Value );
    
    [ misfit1 ] = b4_CALC_MISFIT( TD.Value,predata,par,0,SWwt ); % misfit has structures of summed errors

%% =======================  CALCULATE LIKELIHOOD  =========================
    [ log_likelihood1,misfit1 ] = b5_CALC_LIKELIHOOD( misfit1,TD.Value,model1.datahparm,par);
    
%     fprintf('MISFITS: Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.SpRF,misfit.PsRF,misfit.SW)
%     fprintf('CHI2S:   Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.chi2_sp,misfit.chi2_ps,misfit.chi2_SW)
    
    fail_chain = 0;
    predat_save1 = predata0;

    end % while ifpass
    
%% ========================  ACCEPTANCE CRITERION  ========================
    [ ifaccept ] = b6_IFACCEPT( log_likelihood1,log_likelihood,temp,p_bd*ifpass);        
    
    % ======== PLOT ========  if accept
    if ifaccept && par.inv.verbose && fail_chain==0
        plot_TRUvsPRE( TD.Value,predata);
        if strcmp(projname,'SYNTHETICS')
            plot_MOD_TRUEvsTRIAL( TRUEmodel, model1 );
        end
    end
    
%% ========================  IF ACCEPT ==> CHANGE TO NEW MODEL  =========================
    if ifaccept 
        if par.inv.verbose
            fprintf('  *********************\n  Accepting model! logL:  %.4e ==>  %.4e\n  *********************\n',...
                log_likelihood,log_likelihood1)
        end
        % save new model!
        model = model1;
        log_likelihood = log_likelihood1;
        misfit = misfit1;
        predat_save = predat_save1;
        %% >>> testing
%         n_since_acc = 0;
        %% <<< testing
        
    %% UPDATE KERNEL if needed 
        if newK==true
            [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,par,0,HVK_new);                        
            nchain = 0;
        end
    
    else
        if par.inv.verbose, fprintf('  --FAIL--\n'); end
        if newK, delete_mineos_files(ID,'R'); end
        if newK, delete_mineos_files(ID,'L'); end
        
        %% >>> testing
%         n_since_acc = n_since_acc+1;
        %% <<< testing
    end
    
    % restart-chain if immediate failure
    if isinf(log_likelihood), fail_chain=100; break; end 

    %% >>> testing
%     if n_since_acc > 100
%         1
%     end
    %% <<< testing
    
    
%% =========  reset kernel at end of burn in or after too many iter =======
    resetK = false;
    if newK == false && ii == par.inv.burnin
            fprintf('\n RECALCULATING %s KERNEL - end of burn in\n',chainstr);
            resetK = true;
    end
    if newK == false && nchain > par.inv.maxnkchain
            fprintf('\n RECALCULATING %s KERNEL - chain too long\n',chainstr);
            resetK = true;
    end
    if resetK
        % reset the kernels using the current model
        [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,par,1);
        % need to also reset likelihood and misfit to the new, precise data
        % (likelihood may have been artificially high due to kernel forward
        % calc. approximation - if so, need to undo this, or chain will get
        % stuck once we reset kernels).
        [log_likelihood,misfit] = b8_LIKELIHOOD_RESET(par,predata,TD.Value,Kbase,model.datahparm);
        nchain = 0;
    end


        %% >>> testing
%         accdec(ii) = struct('iter',ii,'logL_current',log_likelihood,'logL_proposed',log_likelihood1,'ifacc',ifaccept,'temp',temp,'preF',p_bd*ifpass,'newk',newK);
        %% <<< testing
 
catch
    if par.inv.verbose, fprintf('  --SOME ERROR--\n'); end
    fail_chain = fail_chain+1;
    continue % skip this model!
end % on try-catch
end % on iterations
%% -------------------------- End iteration  ------------------------------
end % on the fail_chain while...
% ----------
fprintf('\n =========== ENDING ITERATIONS %s ===========\n',chainstr)
% save([resdir,'/chainout/',chainstr],'model0','misfits','allmodels')

misfits_perchain{iii} = misfits;
allmodels_perchain{iii} = allmodels;
if isfield(TD.Value,'SW_Ray')
    SWs_perchain{iii} = preSW;	
end


end % parfor loop
%% ========================================================================
%% ========================================================================
%% ----------------------- End loop on chains  ----------------------------
%% ========================================================================
%% ========================================================================
save([resdir,'/matlab'])
fprintf('Duration of entire run: %.0f s\n',(now-t)*86400)

%% Process results
fprintf('  > Processing results\n')
[misfits_perchain,allmodels_perchain,goodchains] = c1_PROCESS_RESULTS( misfits_perchain,allmodels_perchain,par,1,[resdir,'/modMisfits']);
misfits_perchain_original = misfits_perchain;    % misfits_perchain = misfits_perchain_original;
allmodels_perchain_original = allmodels_perchain;% allmodels_perchain = allmodels_perchain_original;
misfits_perchain = misfits_perchain(goodchains);
allmodels_perchain = allmodels_perchain(goodchains);

[ allmodels_collated ] = collate_allmodels_perchain( allmodels_perchain,par );

[ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels_perchain,[resdir,'/hyperparmtrend.pdf'] );
plot_KNOT_TRENDS( allmodels_perchain,par,[resdir,'/knottrends']  )

posterior = c2_BUILD_POSTERIOR(allmodels_collated,par,zatdep);

fprintf('  > Plotting posterior\n')
plot_MODEL_SUMMARY(posterior,1,[resdir,'/posterior.pdf'])

fprintf('  > Plotting prior vs. posterior\n')
plot_PRIORvsPOSTERIOR(prior,posterior,par,1,[resdir,'/prior2posterior.pdf'])
plot_P2P_recovery(prior,posterior,TRUEmodel,par,1,[resdir,'/mparm_recovery_p2p.pdf'])

fprintf('  > Plotting model suite\n')
[ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels_collated,par );
plot_SUITE_of_MODELS( suite_of_models,posterior,1,[resdir,'/suite_of_models.pdf'],[stadeets(1).Latitude,stadeets(1).Longitude]);
plot_HEATMAP_ALLMODELS(suite_of_models,par,1,[resdir,'/heatmap_of_models.pdf']);

%% Save some things
fprintf('  > Saving misfits, allmods, posterior, model suite\n')
save([resdir,'/misfits_perchain'],'misfits_perchain');
save([resdir,'/allmodels_perchain'],'allmodels_perchain');
save([resdir,'/posterior'],'posterior');
% save([resdir,'/allmods_collated'],'allmodels_collated');
save([resdir,'/mod_suite'],'suite_of_models');
save([resdir,'/goodchains'],'goodchains');
save([resdir,'/SWs_pred'],'SWs_perchain');

%% Final interpolated model with errors
final_model = c4_FINAL_MODEL(posterior,allmodels_collated,par,1,[resdir,'/final_model']);
plot_FINAL_MODEL( final_model,posterior,1,[resdir,'/final_model.pdf'],true,[stadeets(1).Latitude,stadeets(1).Longitude]);

%% predict data with the final model, and calculate the error!
[ final_predata ] = c5_FINAL_FORWARD_MODEL( final_model,par,trudata );

% distribute data for different processing (e.g. _lo, _cms)
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt}; 
    pdt = parse_dtype( dtype ); 
    if strcmp(pdt{1},'BW') && (~strcmp(pdt{3},'def') || ~strcmp(pdt{4},'def'))
        if any(strcmp(par.inv.datatypes,['BW_',pdt{2}])) % only if there IS a standard!
            disp(['replacing ',dtype,' with ',[pdt{1},'_',pdt{2}]])
            final_predata.(dtype) = final_predata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
        end
    end
end
% window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ final_predata ] = predat_process( final_predata,dtype,par);
end
   
[ final_misfit ] = b4_CALC_MISFIT( trudata,final_predata,par,0 );
[ final_log_likelihood,final_misfit ] = b5_CALC_LIKELIHOOD( final_misfit,trudata,final_model.hyperparms,par );
plot_TRUvsPRE( trudata,final_predata,1,[resdir,'/final_true_vs_pred_data.pdf']);
plot_TRUvsPRE_WAVEFORMS( trudata,final_predata,1,[resdir,'/final_true_vs_pred_data_wavs.pdf']);

%addpath([resdir]); 

plot_FIG2_FIT_MODEL( final_model,posterior,prior,par,1,[resdir,'/fig2_FIT_MODEL.pdf']);

% did we save the data?
if ifsavedat
    savedat.gdmods = find([allmodels.bestmods]');
    savedat.gdmods(savedat.gdmods==0) = [];
    save([resdir,'/savedat'],'savedat');
% plot_DATAFITS(trudata,savedat,par,1)
    plot_FIG1_FIT_DATA( trudata,savedat,par,1,[resdir,'/fig1_FIT_DATA.pdf'])
end


%rmpath([resdir]);

% clear('TRUEmodel')
profile viewer

return

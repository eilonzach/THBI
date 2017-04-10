clear
close all

projname = 'WYOMING';
sta = 'RSSD'; % no Ps!!
nwk = 'IU';
gc = 70; % will search for gcarcs +/-3 of this value;
% baz = 315;
global projdir THBIpath
THBIpath = '/Users/zeilon/Documents/MATLAB/BayesianJointInv';
projdir = [THBIpath,'/',projname,'/'];


%% ------------------------- START ------------------------- 
run([THBIpath,'/a0_STARTUP_BAYES']);
cd(projdir);

STAMP=[sta,datestr(now,'_yyyymmddHHMM_pll')];
mkdir([projdir,STAMP]);

%% PARMS
run parms/bayes_inv_parms
par_ORIG = par;
% par came from above script
save([projdir,STAMP,'/par'],'par');
par.inv.verbose=false;

%% PRIOR
fprintf('  > Building prior distribution from %.0f runs\n',par.inv.niter)
prior = a2_BUILD_PRIOR(par,par.inv.niter);
plot_MODEL_SUMMARY(prior,1,[projdir,STAMP,'/prior_fig.pdf']);
save([projdir,STAMP,'/prior'],'prior');
% clone_figure(90,89);

%% stations

%% data info
global TRUEmodel

%% Load & prep data
fprintf('LOADING data\n')
if strcmp(projname,'SYNTHETICS')
    fprintf('Creating custom model and synthetic data\n')
    % z0_SYNTH_MODEL_custommod(par,1); global TRUEmodel; close(95)
    z0_SYNTH_MODEL_simplemod(par,1);  %close(95)
    [ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format

else
    addpath('matguts')
    [~,~,~,TRUEmodel.Z,TRUEmodel.vs,TRUEmodel.vp,TRUEmodel.rho] = RD_1D_Vprofile; close(gcf);
    trudata = load_data(projname,sta,nwk,gc);
    save([projdir,STAMP,'/trudata_ORIG'],'trudata');
end

trudata_ORIG = trudata; trudata_ORIG.PsRF_lo=trudata_ORIG.PsRF; trudata_ORIG.SpRF_lo=trudata_ORIG.SpRF;

% window, filter data 
trudata.PsRF_lo = trudata_ORIG.PsRF;
trudata.SpRF_lo = trudata_ORIG.SpRF;
[ trudata ] = predat_process( trudata,'PsRF',par);
[ trudata ] = predat_process( trudata,'SpRF',par);
[ trudata ] = predat_process( trudata,'PsRF_lo',par);
[ trudata ] = predat_process( trudata,'SpRF_lo',par);
plot_TRUvsPRE_WAVEFORMS(trudata,trudata_ORIG);
plot_TRUvsPRE(trudata,trudata);
    
    
% calc noise parms

%% ---------------------------- INITIATE ----------------------------
profile clear
profile on

%% START DIFFERENT MARKOV CHAINS IN PARALLEL
% nchains = ;
% ifpass = zeros(nchains,1);
model0_perchain = cell(par.inv.nchains,1);
misfits_perchain = cell(par.inv.nchains,1);
allmodels_perchain = cell(par.inv.nchains,1);

%% ========================================================================
%% ========================================================================
fprintf('\n =========== STARTING PARALLEL CHAINS ===========\n')
%% ========================================================================
%% ========================================================================
t = now;
parfor iii = 1:par.inv.nchains 

%% Fail-safe to restart chain if there's a succession of failures
fail_chain=10;
while fail_chain>=10


%% Prep posterior structure
[ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par);

%% Initiate model
ifpass = 0;
% only let starting model out of the loop if it passes conditions
while ifpass==0
    while ifpass==0 % first make sure the starting model satisfies conditions
        model0 = b1_INITIATE_MODEL(par);
        model = model0;
        ifpass = a1_TEST_CONDITIONS( model, par );
    end

    %% starting model kernel
    fprintf('\nCreating starting kernel %s\n',char(64+iii))
    phV = run_mineos(model,trudata.SW.periods,['start',char(64+iii)],0,0,1);
    K  = run_kernels(model,trudata.SW.periods,['start',char(64+iii)],1,0,1);
%         [phV,~,K] = run_mineos_old(model,trudata.SW.periods,[],[],[],1);
    Kbase = struct('modelk',model,'phV',phV,'K',{K});
    
    if any(isnan(phV)),ifpass = false; end
end % now we have a starting model!


%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
ptb = cell({});
nchain = 0;
fail_chain = 0;
ifaccept=true; 
% reset_likelihood;
misfits.lastlogL = -Inf; 
predata=[];
% not parfor
fprintf('\n =========== STARTING ITERATIONS %s ===========\n',char(64+iii))
for ii = 1:par.inv.niter
try
    if rem(ii,20)==0, fprintf('Iteration %s%.0f\n',char(64+iii),ii); end
    ifaccept=false;
    newK = false;
    if fail_chain>9, break, end
    
    % temperature - for perturbation scaling and likelihood increase
    temp = (par.inv.tempmax-1)*erfc(2*(ii-1)./par.inv.cooloff) + 1;
%     if round_level(temp,0.01)>1
%         if par.inv.verbose, fprintf('TEMPERATURE = %.2f\n',temp); end
%     end
    
%% ===========================  PERTURB  ===========================  
    if ii==1
        model1 = model; % don't perturb on first run
        ptbnorm = 0;
        p_bd = 1;
        ptb{ii,1} = 'start';
    else
        ifpass = 0;
        % only let perturbed model out of the loop if it passes conditions
        while ifpass==0
            [model1,ptb{ii,1},p_bd] = b2_PERTURB_MODEL(model,par,temp);
            ifpass = a1_TEST_CONDITIONS( model1, par );
            if ~ifpass, if par.inv.verbose, fprintf('  nope\n'); end; end
            
            [ modptb ] = calc_Vperturbation( Kbase.modelk,model1);
            ptbnorm = norm(modptb.dvsv) + norm(modptb.dvsh);
        end
    end

    if ptbnorm > par.inv.kerneltolmin
        nchain = nchain + 1;
    else
        nchain = nchain + 0.2;
    end

%% ===========================  FORWARD MODEL  ===========================
    % make random run ID (to avoid overwrites in parfor)
	if ~strcmp('sigma',ptb{ii}(end-4:end)) || isempty(predata)
		id = [char(64+iii),num2str(round(1e9*(now-t))),num2str(randi(9)),num2str(randi(9))];
	
		predata = b3_FORWARD_MODEL( model1,Kbase,par,trudata,id,0); predata0=predata;
    
        % continue if any Sp or PS inhomogeneous or weird output
        if predata.PsRF.nsamp<predata.PsRF.samprate*diff(par.datprocess.Twin.PsRF)
            fprintf('Not enough P data!\n'),fail_chain=fail_chain+1;continue, end
        if predata.SpRF.nsamp<predata.SpRF.samprate*diff(par.datprocess.Twin.SpRF)
            fprintf('Not enough S data!\n'),fail_chain=fail_chain+1;continue, end
        if any(any(isnan(predata.SpRF.ZRT))) || any(any(isnan(predata.PsRF.ZRT))) 
            fprintf('inhomogeneous!\n'), fail_chain=fail_chain+1; continue, end
        
        for idt = 1:length(par.inv.datatypes)
            [ predata ] = predat_process( predata,par.inv.datatypes{idt},par);
        end
        
        % reject if a nan channel
        if any(any(isnan(predata.PsRF.ZRT))) ||  any(any(isnan(predata.SpRF.ZRT))), 
            fprintf('NaN DATA!\n'),fail_chain=fail_chain+1;continue, end
	
		% Explicitly use mineos if ptb is too large
		if par.inv.verbose, fprintf('    Perturbation %.2f\n',ptbnorm); end
		if ptbnorm/par.inv.kerneltolmax > random('unif',par.inv.kerneltolmin/par.inv.kerneltolmax,1,1) %#ok<ALIGN> % control chance of going to MINEOS
			swk = predata.SW.phV; % record existing phV from kernels
		
			predata.SW.phV = run_mineos(model1,trudata.SW.periods,id,0,0,par.inv.verbose);
		
			if par.inv.verbose
				fprintf('    RMS diff is %.4f\n',rms(swk-predata.SW.phV)); % RMS difference
			end
		
			newK = true;
		else 
			Ktry = [];
		end
	end % only redo data if model has changed 

%      plot_TRUvsPRE_old(trudata,predata)]

        % continue if any Sp or PS inhomogeneous or weird output
        if predata.PsRF.nsamp<predata.PsRF.samprate*diff(par.datprocess.Twin.PsRF)
            fprintf('Not enough P data!\n'),fail_chain=fail_chain+1;continue, end
        if predata.SpRF.nsamp<predata.SpRF.samprate*diff(par.datprocess.Twin.SpRF)
            fprintf('Not enough S data!\n'),fail_chain=fail_chain+1;continue, end
        if any(any(isnan(predata.SpRF.ZRT))) || any(any(isnan(predata.PsRF.ZRT))) 
            fprintf('inhomogeneous!\n'), fail_chain=fail_chain+1; continue, end
    
%% =========================  CALCULATE MISFIT  ===========================
    
    % if weighting SWs:
    if all(par.inv.Kweight == true) % if using default weight by fraction of kernel in model
        SWwt = calc_K_in_model( Kbase.K,par );
        SWwt=SWwt/mean(SWwt);
    elseif all(par.inv.Kweight == false) % if no weight
        SWwt = [];
    else
        SWwt = par.inv.Kweight; % if not explicitly "false", assume custom wts and use.
    end
    
    [ misfit ] = b4_CALC_MISFIT( trudata,predata,par,0,SWwt ); % misfit has structures of summed errors

%% =======================  CALCULATE LIKELIHOOD  =========================
    [ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par);
    
%     fprintf('MISFITS: Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.SpRF,misfit.PsRF,misfit.SW)
%     fprintf('CHI2S:   Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.chi2_sp,misfit.chi2_ps,misfit.chi2_SW)
    
%% ========================  ACCEPTANCE CRITERION  ========================
    [ ifaccept ] = b6_IFACCEPT( log_likelihood+log(p_bd),misfits,temp );
    
    % ======== PLOT ========  if accept
    if ifaccept && par.inv.verbose
        plot_TRUvsPRE( trudata,predata);
        if strcmp(projname,'SYNTHETICS')
            plot_MOD_TRUEvsTRIAL( TRUEmodel, model1 );
        end
    end
    
%% ========================  IF ACCEPT ==> STORE  =========================
    if ifaccept 
        if par.inv.verbose
            fprintf('  *********************\n  Accepting model! logL:  %.4e ==>  %.4e\n  *********************\n',...
                misfits.lastlogL,log_likelihood)
        end
        model = model1;
        
        [misfits,allmodels,savedat] = b9_SAVE_RESULT(ii,log_likelihood,misfit,model1,misfits,allmodels,predata0,savedat);
    %% UPDATE KERNEL if needed 
        if newK==true
            Kbase = struct('modelk',model);
            Kbase.phV = predata.SW.phV;
            Kbase.K = run_kernels(model,trudata.SW.periods,id,1,0,par.inv.verbose);
            nchain = 0;
        end
    
    else
        if par.inv.verbose, fprintf('  --FAIL--\n'); end
        if newK, delete_mineos_files(id,trudata.SW.periods); end
    end
    
    fail_chain = 0;
    
%% =========  redo kernel at end of burn in or if chain is too long =======
    if (newK == false) && (nchain > 30);
        if par.inv.verbose, fprintf('\n RECALCULATING KERNEL - too long chain\n'); end
        Kbase.modelk = model;
        Kbase.phV = run_mineos(model,trudata.SW.periods,id,0,0,par.inv.verbose);
        Kbase.K  = run_kernels(model,trudata.SW.periods,id,1,0,par.inv.verbose);
        nchain = 0;
    end
    if ii == par.inv.burnin
%         fprintf('\n RECALCULATING KERNEL - end of burn in\n')
        Kbase.modelk = model;
        Kbase.phV = run_mineos(model,trudata.SW.periods,id,0,0,par.inv.verbose);
        Kbase.K  = run_kernels(model,trudata.SW.periods,id,1,0,par.inv.verbose);
        nchain = 0;
    end
    
catch
    if par.inv.verbose, fprintf('  --SOME ERROR--\n'); end
    fail_chain = fail_chain+1;
    continue % skip this model!
end % on try-catch
end % oniterations
%% -------------------------- End iteration  ------------------------------
end % on the fail_chain while...
% ----------

model0_perchain{iii} = model0;
misfits_perchain{iii} = misfits;
allmodels_perchain{iii} = allmodels;


end % parfor loop
%% ========================================================================
%% ========================================================================
%% ----------------------- End loop on chains  ----------------------------
%% ========================================================================
%% ========================================================================
save([STAMP,'/matlab'])
fprintf('Duration of entire run: %.0f s\n',(now-t)*86400)

%% Process results
fprintf('  > Processing results\n')
[misfits_perchain,allmodels_perchain,goodchains] = c1_PROCESS_RESULTS( misfits_perchain,allmodels_perchain,par,1,['figs/misfits_vs_iter_',STAMP]);
misfits_perchain_original = misfits_perchain;
allmodels_perchain_original = allmodels_perchain;
misfits_perchain = misfits_perchain(goodchains);
allmodels_perchain = allmodels_perchain(goodchains);

[ allmodels_collated ] = collate_allmodels_perchain( allmodels_perchain,par );
[ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels_perchain );
plot_KNOT_TRENDS( allmodels_perchain,par,[projdir,STAMP,'/knottrends.pdf']  )

posterior = c2_BUILD_POSTERIOR(allmodels_collated,par);
fprintf('  > Saving posterior\n')
save([projdir,STAMP,'/posterior'],'posterior');
save([projdir,STAMP,'/allmods_collated'],'allmodels_collated');

fprintf('  > Plotting posterior\n')
plot_MODEL_SUMMARY(posterior,1,[projdir,STAMP,'/posteror.pdf'])

fprintf('  > Plotting prior vs. posterior\n')
load([projdir,STAMP,'/prior']) ; load([projdir,STAMP,'/posterior'])
plot_PRIORvsPOSTERIOR(prior,posterior,par,1,[projdir,STAMP,'/prior2posterior.pdf'])

[ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels_collated,par );
fprintf('  > Saving model suite\n')
save([projdir,STAMP,'/mod_suite'],'suite_of_models');
fprintf('  > Plotting model suite\n')
plot_SUITE_of_MODELS( suite_of_models,posterior,1,[projdir,STAMP,'/suite_of_models.pdf'])


%% Final interpolated model with errors
final_model = c4_FINAL_MODEL(posterior,allmodels_collated,par,1,[projdir,STAMP,'/final_model']);
plot_FINAL_MODEL( final_model,posterior,1,[projdir,STAMP,'/final_model.pdf']);

%% predict data with the final model, and calculate the error!
[ final_predata ] = final_forward_model( final_model,par,trudata );
final_predata.PsRF_lo = final_predata.PsRF;
final_predata.SpRF_lo = final_predata.SpRF;
[ final_predata ] = predat_process( final_predata,'PsRF',par);
[ final_predata ] = predat_process( final_predata,'SpRF',par);
[ final_predata ] = predat_process( final_predata,'PsRF_lo',par);
[ final_predata ] = predat_process( final_predata,'SpRF_lo',par);

[ final_misfit ] = b4_CALC_MISFIT( trudata,final_predata,par,0 );
[ final_log_likelihood,final_misfit ] = b5_CALC_LIKELIHOOD( final_misfit,trudata,...
                 struct('sigmaSpRF',10^final_model.hyperparms.sigmaSpRF.mu_log10,...
                        'sigmaPsRF',10^final_model.hyperparms.sigmaPsRF.mu_log10,...
                        'sigmaSpRF_lo',10^final_model.hyperparms.sigmaSpRF_lo.mu_log10,...
                        'sigmaPsRF_lo',10^final_model.hyperparms.sigmaPsRF_lo.mu_log10,...
                        'sigmaSW',10^final_model.hyperparms.sigmaSW.mu_log10),par );
plot_TRUvsPRE( trudata,final_predata,1,[projdir,STAMP,'/final_true_vs_pred_data_fig.pdf']);

addpath([projdir,STAMP]); 
plot_FIG2_FIT_MODEL(final_model,posterior,prior,par,1,[projdir,STAMP,'/fig2_MODEL.pdf']);
rmpath([projdir,STAMP]);

clear('TRUEmodel')
profile viewer

return

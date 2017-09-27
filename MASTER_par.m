clear all
close all


projname = 'WYOMING'; % SYNTHETICS or WYOMING, for now
sta = 'EYMN';
nwk = 'US';
gc = [73,69,66,59,44,38,32]; % will search for gcarcs +/-3 of this value;
datN = 20;
% baz = 315;

notes = [...
    'Par inversion with data made using new clustering scheme \n',...
    'bining by baz and gcarc.  \n'...
    'Using all data types: Ps,Pslo,Pscms,Sp,Splo,SW.\n',...
    'This run has a sediment layer DISallowed.\n',...
    'Max depth 300 km.\n',...
        ];

%% ------------------------- START ------------------------- 
global projdir THBIpath TRUEmodel
THBIpath = '/Users/zeilon/Documents/MATLAB/BayesianJointInv';
projdir = [THBIpath,'/',projname,'/'];
avardir = sprintf('STA_inversions/%s_dat%.0f/',sta,datN);
stadeets = irisFetch.Stations('station',nwk,sta,'*','*');

run([THBIpath,'/a0_STARTUP_BAYES']);
cd(projdir);

%% PARMS
run parms/bayes_inv_parms
if strcmp(projname,'SYNTHETICS'),
    if strcmp(par.synth.noisetype,'real'), sta=['SYNTH_',sta]; else sta = 'SYNTH'; end
    par.sta = sta; 
    par.nwk = '--';
else
    par.sta = sta;
    par.nwk = nwk;
    par.gc = gc;
end
par.inv.verbose=false;

STAMP=[sta,datestr(now,'_yyyymmddHHMM_pll')]; par.STAMP = STAMP;
resdir = [projdir,avardir,STAMP];
mkdir(resdir);
fid = fopen([resdir,'/notes.txt'],'w'); fprintf(fid,notes); fclose(fid);

par_ORIG = par;
% par came from above script
save([resdir,'/par'],'par');
copyfile('parms/bayes_inv_parms.m',resdir);

for id = 1:length(par.inv.datatypes);
    allpdytp(id,:)=parse_dtype(par.inv.datatypes{id});
end

%% Load & prep data
fprintf('LOADING data\n')
if strcmp(projname,'WYOMING')

    addpath('matguts')
%     [~,~,~,TRUEmodel.Z,TRUEmodel.vs,TRUEmodel.vp,TRUEmodel.rho] = RD_1D_Vprofile; close(gcf);
    [trudata,cheatstr] = load_data(avardir,sta,nwk,gc);
    sta = [sta,cheatstr];
    % distribute data for different processing (e.g. _lo, _cms)
    for idt = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{idt}; pdt = parse_dtype( dtype ); 
        if ~isfield(trudata,par.inv.datatypes{idt}) && strcmp(pdt{1},'BW') 
            trudata.(dtype) = trudata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
        end
    end
    
elseif strcmp(projname,'SYNTHETICS')
    
    fprintf(' > Creating custom model and synthetic data\n')
    
    % make synth model
    addpath('~/Documents/MATLAB/THBI_paper/Figure_3/')
    z0_SYNTH_MODEL_fig3_splinemod(par,0);  %close(95)
    save([resdir,'/trumodel'],'TRUEmodel');

    % make synth data
    [ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
    if strcmp(par.synth.noisetype,'real')
        noise_sta_deets = struct('datadir',[THBIpath,'/WYOMING/AVARS/'],'sta',sta,'nwk',nwk,'gc',par.synth.gcarcs);
        par.synth.noise_sta_deets = noise_sta_deets;
        [ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,noise_sta_deets );
    end 
end

% save data
save([resdir,'/trudata_ORIG'],'trudata');

% window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ trudata ] = predat_process( trudata,dtype,par);
end
plot_TRU_WAVEFORMS(trudata);
% plot_TRUvsPRE_WAVEFORMS(trudata,trudata_ORIG);
plot_TRUvsPRE(trudata,trudata);

%% PRIOR
fprintf('  > Building prior distribution from %.0f runs\n',min([par.inv.niter,1e4]))
prior = a2_BUILD_PRIOR(par,min([par.inv.niter,1e4]));
plot_MODEL_SUMMARY(prior,1,[resdir,'/prior_fig.pdf']);
save([resdir,'/prior'],'prior');
    
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
        model = expand_dathparms_to_data( model,trudata,par );
        ifpass = a1_TEST_CONDITIONS( model, par );
    end

    %% starting model kernel
    fprintf('\nCreating starting kernel %s\n',char(64+iii))

    Kbase = initiate_Kbase; Kbase.modelk = model;
    for id = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{id};
        pdtyp = parse_dtype(dtype); 
        if ~strcmp(pdtyp{1},'SW'), continue; end
        [phV,grV] = run_mineos(model,trudata.(dtype).periods,pdtyp{2},['start',char(64+iii)],0,0,1);
        K  = run_kernels(model,trudata.(dtype).periods,pdtyp{2},pdtyp{3},['start',char(64+iii)],1,0,1);
        [ Kbase ] = populate_Kbase( Kbase,dtype,phV,grV,{K} );
    
        if any(isnan(phV)),ifpass = false; end
    end % loop on datatypes

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
    
    nchain  = kchain_addcount( nchain,ptbnorm,par );
    

%% ===========================  FORWARD MODEL  ===========================
    % make random run ID (to avoid overwrites in parfor)
	if ~strcmp('sigma',ptb{ii}(end-4:end)) || isempty(predata)
		ID = [char(64+iii),num2str(round(1e9*(now-t))),num2str(randi(9)),num2str(randi(9))];

		try
            predata = b3_FORWARD_MODEL( model1,Kbase,par,trudata,ID,0); predata0=predata;
        catch
            fprintf('Forward model error\n'); fail_chain=fail_chain+1; continue
        end
        
        % continue if any Sp or PS inhomogeneous or nan or weird output
        if ifforwardfail(predata,par), fail_chain=fail_chain+1; continue, end
        
        for idt = 1:length(par.inv.datatypes)
            [ predata ] = predat_process( predata,par.inv.datatypes{idt},par);
        end
        
		% Explicitly use mineos if ptb is too large
		if par.inv.verbose, fprintf('    Perturbation %.2f\n',ptbnorm); end
		if ptbnorm/par.inv.kerneltolmax > random('unif',par.inv.kerneltolmed/par.inv.kerneltolmax,1,1) % control chance of going to MINEOS
            SW = [];
            if any(strcmp(allpdytp(:,2),'Ray')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'Ray'),1,'first'));
                [SW.Ray.phV,SW.Ray.grV] = run_mineos(model1,trudata.(itp{1}).periods,'R',ID,0,0,par.inv.verbose);
            end
            if any(strcmp(allpdytp(:,2),'Lov')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'Lov'),1,'first'));
                [SW.Lov.phV,SW.Lov.grV] = run_mineos(model1,trudata.(itp{1}).periods,'L',ID,0,0,par.inv.verbose);
            end
            for id = 1:length(par.inv.datatypes)
                dtype = par.inv.datatypes{id}; pdtyp=parse_dtype(dtype); if ~strcmp(pdtyp{1},'SW'),continue; end
                swk = predata.(dtype).phV; % record existing phV from kernels

                predata.(dtype).phV = SW.(pdtyp{2}).(pdtyp{3});

                if par.inv.verbose
                    fprintf('    RMS diff is %.4f\n',rms(swk-predata.(dtype).phV)); % RMS difference
                end
            end
			newK = true;
		else 
			Ktry = [];
		end
    end % only redo data if model has changed 

%      plot_TRUvsPRE_old(trudata,predata)]

    % continue if any Sp or PS inhomogeneous or nan or weird output
    if ifforwardfail(predata,par), fail_chain=fail_chain+1; continue, end

    
%% =========================  CALCULATE MISFIT  ===========================
    
    % SW weights, if applicable 
    [ SWwt ] = make_SW_weight( par,Kbase );
    
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
            Kbase.modelk = model;
            for id = 1:length(par.inv.datatypes)
                dtype = par.inv.datatypes{id};pdtyp = parse_dtype(dtype); if ~strcmp(pdtyp{1},'SW'), continue; end
                K = run_kernels(model,trudata.(dtype).periods,pdtyp{2},pdtyp{3},ID,1,0,par.inv.verbose);
                Kbase = populate_Kbase( Kbase,dtype,predata.(dtype).phV,[],{K} );
            end
            nchain = 0;
        end
    
    else
        if par.inv.verbose, fprintf('  --FAIL--\n'); end
        if newK, delete_mineos_files(ID,'R'); end
        if newK, delete_mineos_files(ID,'L'); end
    end
    
    fail_chain = 0;
    
%% =========  redo kernel at end of burn in or if chain is too long =======
    if (newK == false) && (nchain > par.inv.maxnkchain);
        if par.inv.verbose, fprintf('\n RECALCULATING KERNEL - too long chain\n'); end
        Kbase.modelk = model;
        for id = 1:length(par.inv.datatypes)
            dtype = par.inv.datatypes{id};pdtyp = parse_dtype(dtype); if ~strcmp(pdtyp{1},'SW'), continue; end
            [phV,grV] = run_mineos(model,trudata.(dtype).periods,pdtyp{2},ID,0,0,par.inv.verbose);
            K = run_kernels(model,trudata.(dtype).periods,pdtyp{2},pdtyp{3},ID,1,0,par.inv.verbose);
            Kbase = populate_Kbase( Kbase,dtype,phV,grV,{K} );
        end
        nchain = 0;
    end
    if ii == par.inv.burnin
%         fprintf('\n RECALCULATING KERNEL - end of burn in\n')
        Kbase.modelk = model;
        for id = 1:length(par.inv.datatypes)
            dtype = par.inv.datatypes{id};pdtyp = parse_dtype(dtype); if ~strcmp(pdtyp{1},'SW'), continue; end
            [phV,grV] = run_mineos(model,trudata.(dtype).periods,pdtyp{2},ID,0,0,par.inv.verbose);
            K = run_kernels(model,trudata.(dtype).periods,pdtyp{2},pdtyp{3},ID,1,0,par.inv.verbose);
            Kbase = populate_Kbase( Kbase,dtype,phV,grV,{K} );
        end
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
fprintf('\n =========== ENDING ITERATIONS %s ===========\n',char(64+iii))

model0_perchain{iii} = model0;
misfits_perchain{iii} = misfits;
allmodels_perchain{iii} = allmodels;


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
misfits_perchain_original = misfits_perchain;
allmodels_perchain_original = allmodels_perchain;
misfits_perchain = misfits_perchain(goodchains);
allmodels_perchain = allmodels_perchain(goodchains);

[ allmodels_collated ] = collate_allmodels_perchain( allmodels_perchain,par );

[ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels_perchain,[resdir,'/hyperparmtrend.pdf'] );
plot_KNOT_TRENDS( allmodels_perchain,par,[resdir,'/knottrends']  )

posterior = c2_BUILD_POSTERIOR(allmodels_collated,par);

fprintf('  > Plotting posterior\n')
plot_MODEL_SUMMARY(posterior,1,[resdir,'/posterior.pdf'])

fprintf('  > Plotting prior vs. posterior\n')
plot_PRIORvsPOSTERIOR(prior,posterior,par,1,[resdir,'/prior2posterior.pdf'])

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

%% Final interpolated model with errors
final_model = c4_FINAL_MODEL(posterior,allmodels_collated,par,1,[resdir,'/final_model']);
plot_FINAL_MODEL( final_model,posterior,1,[resdir,'/final_model.pdf'],true,[stadeets(1).Latitude,stadeets(1).Longitude]);

%% predict data with the final model, and calculate the error!
[ final_predata ] = c5_FINAL_FORWARD_MODEL( final_model,par,trudata );

% distribute data for different processing (e.g. _lo, _cms)
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt}; pdt = parse_dtype( dtype ); 
    if ~isfield(final_predata,par.inv.datatypes{idt}) && strcmp(pdt{1},'BW') 
        final_predata.(dtype) = final_predata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
    end
end
% window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ final_predata ] = predat_process( final_predata,dtype,par);
end
   
[ final_misfit ] = b4_CALC_MISFIT( trudata,final_predata,par,0 );
% [ final_log_likelihood,final_misfit ] = b5_CALC_LIKELIHOOD( final_misfit,trudata,...
%                  struct('sigmaSpRF',10.^final_model.hyperparms.sigmaSpRF.mu_log10,...
%                         'sigmaPsRF',10.^final_model.hyperparms.sigmaPsRF.mu_log10,...
%                         'sigmaSpRF_lo',10.^final_model.hyperparms.sigmaSpRF_lo.mu_log10,...
%                         'sigmaPsRF_lo',10.^final_model.hyperparms.sigmaPsRF_lo.mu_log10,...
%                         'sigmaSW',10.^final_model.hyperparms.sigmaSW.mu_log10),par );
plot_TRUvsPRE( trudata,final_predata,1,[resdir,'/final_true_vs_pred_data.pdf']);
plot_TRUvsPRE_WAVEFORMS( trudata,final_predata,1,[resdir,'/final_true_vs_pred_data_wavs.pdf']);

%addpath([resdir]); 

plot_FIG2_FIT_MODEL( final_model,posterior,prior,par,1,[resdir,'/fig2_FIT_MODEL.pdf']);

%rmpath([resdir]);

clear('TRUEmodel')
profile viewer

return

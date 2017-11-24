clear all
close all


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
projdir = [THBIpath,'/test_detailed_balance/'];

run([THBIpath,'/a0_STARTUP_BAYES']);
cd(projdir);

%% PARMS
run('bayes_inv_parms')

par.inv.verbose=false;
par.sta = 'detailed_balance'; par.nwk=[];

STAMP=['DetBal',datestr(now,'_yyyymmddHHMM_pll')]; par.STAMP = STAMP;
resdir = projdir;
mkdir(resdir);
fid = fopen([resdir,'/notes.txt'],'w'); fprintf(fid,notes); fclose(fid);

par_ORIG = par;
% par came from above script
save([resdir,'/par_',STAMP],'par');
copyfile('bayes_inv_parms.m',[resdir,STAMP,'_parms.m']);

for id = 1:length(par.inv.datatypes)
    allpdytp(id,:)=parse_dtype(par.inv.datatypes{id});
end

%% PRIOR
fprintf('  > Building prior distribution from %.0f runs\n',min([par.inv.niter,1e4]))
prior = a2_BUILD_PRIOR(par,min([par.inv.niter,1e4]));
plot_MODEL_SUMMARY(prior,1,[resdir,'/prior_fig.pdf']);
save([resdir,'/prior'],'prior');
    
%% ---------------------------- INITIATE ----------------------------
profile clear
profile on

%% START DIFFERENT MARKOV CHAINS IN PARALLEL
model0_perchain = cell(par.inv.nchains,1);
allmodels_perchain = cell(par.inv.nchains,1);

%% ========================================================================
%% ========================================================================
fprintf('\n =========== STARTING PARALLEL CHAINS ===========\n')
%% ========================================================================
%% ========================================================================
t = now;
parfor iii = 1:par.inv.nchains 

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
end % now we have a starting model!


%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
ptb = cell({});
nchain = 0;
ifaccept=true; 
% reset_likelihood;
misfits.lastlogL = -Inf; 
predata=[];
% not parfor
fprintf('\n =========== STARTING ITERATIONS %s ===========\n',char(64+iii))
for ii = 1:par.inv.niter
    if rem(ii,20)==0, fprintf('Iteration %s%.0f\n',char(64+iii),ii); end
    ifaccept=false;
    newK = false;
    
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
            if ~ifpass, 
                if par.inv.verbose, fprintf('  nope\n'); end; 
            end
        end
    end
    
%% ========================  ACCEPTANCE CRITERION  ========================
    log_likelihood = log(p_bd);
    [ ifaccept ] = b6_IFACCEPT( log_likelihood,misfits,temp );
        
%% ========================  IF ACCEPT ==> STORE  =========================
    if ifaccept 
        if par.inv.verbose
            fprintf('  *********************\n  Accepting model! logL:  %.4e ==>  %.4e\n  *********************\n',...
                misfits.lastlogL,log_likelihood)
        end
        model = model1;
        
        [allmodels] = DB_save_result(ii,log_likelihood,model1,allmodels);

    else
        if par.inv.verbose, fprintf('  --FAIL--\n'); end
    end
    
    
end % on iterations
%% -------------------------- End iteration  ------------------------------
fprintf('\n =========== ENDING ITERATIONS %s ===========\n',char(64+iii))

model0_perchain{iii} = model0;
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
[ allmodels_collated ] = collate_allmodels_perchain( allmodels_perchain,par );
allmodels_collated = dealto(allmodels_collated,'bestmods',true(length(allmodels_collated),1));
for iii = 1:par.inv.nchains, allmodels_perchain{iii} = dealto(allmodels_perchain{iii},'bestmods',true(length(allmodels_perchain{iii}),1)); end
% [ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels_perchain,[resdir,'/hyperparmtrend.pdf'] );
plot_KNOT_TRENDS( allmodels_perchain,par,[resdir,'/knottrends']  )
load prior
posterior = c2_BUILD_POSTERIOR(allmodels_collated,par);

fprintf('  > Plotting posterior\n')
plot_MODEL_SUMMARY(posterior,1,[resdir,'/posterior.pdf'])

fprintf('  > Plotting prior vs. posterior\n')
plot_PRIORvsPOSTERIOR(prior,posterior,par,1,[resdir,'/prior2posterior.pdf'])

fprintf('  > Plotting model suite\n')
[ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels_collated,par );
plot_SUITE_of_MODELS( suite_of_models,posterior,1,[resdir,'/suite_of_models.pdf']);
plot_HEATMAP_ALLMODELS(suite_of_models,par,1,[resdir,'/heatmap_of_models.pdf']);

%% Save some things
fprintf('  > Saving misfits, allmods, posterior, model suite\n')
save([resdir,'/allmodels_perchain'],'allmodels_perchain','-v7.3');
save([resdir,'/posterior'],'posterior');
% save([resdir,'/allmods_collated'],'allmodels_collated');
save([resdir,'/mod_suite'],'suite_of_models');

profile viewer

return

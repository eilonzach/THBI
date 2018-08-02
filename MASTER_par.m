clear all
close all


projname = 'SYNTHETICS'; % SYNTHETICS, LAB_tests, or WYOMING, for now
sta = 'EYMN';
nwk = 'US';
gc = [69,59,40,38,36,66]; % will search for gcarcs +/-3 of this value;
% baz = 315;

notes = [...
    'Par inversion with data made using new clustering scheme \n',...
    'bining by baz and gcarc.  \n'...
    'Using all data types: Ps,Pslo,Pscms,Sp,Splo,SW.\n',...
    'This run has a sediment layer DISallowed.\n',...
    'Max depth 300 km.\n',...
    'New techniques of ignoring inhomog. P conversion layers\n',...
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
	dtps = {'BW_Ps_cms','SW_Ray_phV'};     

	% noise details, if "real"
	noisesta = 'RSSD';
	noisenwk = 'IU';
	noisegcarcs = [73,38];
	noiseshape = 'real'; % 'white' or 'real'
	noiseup = 0.5; % factor to increase real noise

	% naming convention
	dtpstr='_';
	if any(strcmp(dtps,'BW_Ps')), dtpstr=[dtpstr,'Ps']; end
	if any(strcmp(dtps,'BW_Ps_cms')), dtpstr=[dtpstr,'Pscms']; end
	if any(strcmp(dtps,'BW_Sp')), dtpstr=[dtpstr,'Sp']; end
	if any(strcmp(dtps,'SW_Ray_phV')), dtpstr=[dtpstr,'SW']; end

	ifsavedat = true;

	sta = ['LAB_s',num2str(zsed),'_m',num2str(zmoh),'_z',num2str(zlab),'_w',num2str(wlab),'_f',num2str(100*flab),dtpstr];
    par.sta = sta; par.nwk = 'LAB_test';
else
    par.sta = sta;
    par.nwk = nwk;
    par.gc = gc;
    datN = 20;
end

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

%% Load & prep data
fprintf('LOADING data\n')
if strcmp(projname,'SYNTHETICS')
    
    fprintf(' > Creating custom model and synthetic data\n')
    
    % make synth model
    z0_SYNTH_MODEL_fig3_splinemod(par,0);  %close(95)
    save([resdir,'/trumodel'],'TRUEmodel');

    % make synth data
    [ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
    if strcmp(par.synth.noisetype,'real')
        noise_sta_deets = struct('datadir',['/Volumes/data/THBI/US/STAsinv/',noisesta,'_dat20/'],...
								 'sta',noisesta,'nwk',noisenwk,'gc',noisegcarcs,'noiseup',noiseup,'noiseshape',noiseshape);
        par.synth.noise_sta_deets = noise_sta_deets;
        [ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,noise_sta_deets );
    end

    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);

elseif strcmp(projname,'LAB_tests')
	z0_SYNTH_MODEL_LAB_TEST(par,zsed,zmoh,zlab,wlab,flab,1) ;
	save([resdir,'/trumodel'],'TRUEmodel');

	% make synth data
	[ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
	trudata_noiseless = trudata;

	trudata = trudata_noiseless;
	if strcmp(par.synth.noisetype,'real')
		noise_sta_deets = struct('datadir',['/Volumes/data/THBI/US/STAsinv/',noisesta,'_dat20/'],...
								 'sta',noisesta,'nwk',noisenwk,'gc',noisegcarcs,'noiseup',noiseup,'noiseshape',noiseshape);
		par.synth.noise_sta_deets = noise_sta_deets;
	%     [ trudata,par ] = z2_NOISIFY_SYNTH_makestack( trudata, par,noise_sta_deets );
		[ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,noise_sta_deets );
    end
    
    trudata = duplicate_data_distribute(trudata,par);
    stadeets = struct('Latitude',[],'Longitude',[]);

else
	try stadeets = irisFetch.Stations('station',nwk,sta,'*','*'); 
    catch, load([proj.rawdatadir,'/stainfo_master.mat']); 
        stadeets = struct('Latitude',stainfo(strcmp({stainfo.StationCode},sta)).Latitude,...
                          'Longitude',stainfo(strcmp({stainfo.StationCode},sta)).Longitude);
    end

%     [~,~,~,TRUEmodel.Z,TRUEmodel.vs,TRUEmodel.vp,TRUEmodel.rho] = RD_1D_Vprofile; close(gcf);
    [trudata,zeroDstr] = load_data(avardir,sta,nwk,gc);
    sta = [sta,zeroDstr];
    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);
    
    for idt = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{idt}; pdt = parse_dtype( dtype ); 
        if ~isfield(trudata,par.inv.datatypes{idt}) && strcmp(pdt{1},'BW') 
            trudata.(dtype) = trudata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
        end
        % set prior sigma as geometric mean of data sigma
        if isfield(trudata.(dtype),'sigma') && ~isnan(geomean(trudata.(dtype).sigma))
            par.mod.data.prior_sigma.(pdt{1}).(pdt{2}).(pdt{3}) = geomean(trudata.(dtype).sigma);
        end
    end

    % get rid of data that wont be used in inversion - NB NEED EXACT DATA MATCH
    trudtypes = fieldnames(trudata);
    for idt = 1:length(trudtypes)
        if all(~strcmp(trudtypes{idt},par.inv.datatypes)) % no match
            fprintf('WARNING - removing %s data from trudata\n',trudtypes{idt})
            trudata = rmfield(trudata,trudtypes{idt});
        end
        if par.inv.BWclust~=0 && any(regexp(trudtypes{idt},'BW'))
            fprintf('WARNING - removing %s data not in cluster %.0f\n',trudtypes{idt},par.inv.BWclust)
            trudata.(trudtypes{idt}) = trudata.(trudtypes{idt})([trudata.(trudtypes{idt}).clust]==par.inv.BWclust);
        end
    end

end

% save data
save([resdir,'/trudata_ORIG'],'trudata');

% window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ trudata ] = predat_process( trudata,dtype,par);
end
% plot_TRU_WAVEFORMS(trudata);
% plot_TRUvsPRE(trudata,trudata);
save([resdir,'/trudata_USE'],'trudata');


%% PRIOR
fprintf('  > Building prior distribution from %.0f runs\n',max([par.inv.niter,1e5]))
zatdep = [5:5:par.mod.maxz]';
prior = a2_BUILD_EMPIRICAL_PRIOR(par,max([par.inv.niter,1e5]),14,zatdep);
plot_MODEL_SUMMARY(prior,1,[resdir,'/prior_fig.pdf']);
save([resdir,'/prior'],'prior','par');

%% ---------------------------- INITIATE ----------------------------
profile clear
profile on
delete(gcp('nocreate'));

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
        model = expand_dathparms_to_data( model,trudata,par );
        ifpass = a1_TEST_CONDITIONS( model, par );
    end

    %% starting model kernel
    fprintf('\nCreating starting kernels %s\n',chainstr)
    try
        [Kbase] = make_allkernels(model,[],trudata,['start',chainstr],par);
    catch
        ifpass = false; continue; 
    end
    
end % now we have a starting model!


%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
ptb = cell({});
nchain = 0;
fail_chain = 0;
ifaccept=true;
if isfield(trudata,'SW_Ray')
    preSW = zeros(length(trudata.SW_Ray.periods),ceil(par.inv.niter./par.inv.saveperN));
end
% reset_likelihood;
log_likelihood = -Inf;
predata=[]; predat_save = []; misfit = [];
% not parfor
fprintf('\n =========== STARTING ITERATIONS %s ===========\n',chainstr)
for ii = 1:par.inv.niter
    
%% SAVE every saveperN
if mod(ii,par.inv.saveperN)==0 && log_likelihood ~= -Inf
    [misfits,allmodels,savedat] = b9_SAVE_RESULT(ii,log_likelihood,misfit,model,misfits,allmodels,predat_save,savedat);
    if isfield(trudata,'SW_Ray')
        preSW(:,misfits.Nstored) = predata.SW_Ray.phV;
    end
end
    
try
    if rem(ii,par.inv.saveperN)==0 || ii==1, fprintf('Iteration %s%.0f\n',chainstr,ii); end
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
		ptbnorm = norm(modptb.dvsv) + norm(modptb.dvsh);
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
            predata = b3_FORWARD_MODEL( model1,Kbase,par,trudata,ID,0);
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
            [ predata ] = predat_process( predata,par.inv.datatypes{idt},par);
        end
        
		% Explicitly use mineos + Tanimoto scripts if ptb is too large
		if par.inv.verbose, fprintf('    Perturbation %.2f\n',ptbnorm); end
		if ptbnorm/par.inv.kerneltolmax > random('unif',par.inv.kerneltolmed/par.inv.kerneltolmax,1,1) % control chance of going to MINEOS
            SW = struct('Ray',struct('phV',[],'grV',[]),'Lov',struct('phV',[],'grV',[]),'HV',struct('HVr',[]));
            
            if any(strcmp(allpdytp(:,2),'Ray')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'Ray'),1,'first')); %#ok<PFBNS>
                [SW.Ray.phV,SW.Ray.grV] = run_mineos(model1,trudata.(itp{1}).periods,'R',ID,0,0,par.inv.verbose);
            end
            if any(strcmp(allpdytp(:,2),'Lov')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'Lov'),1,'first'));
                [SW.Lov.phV,SW.Lov.grV] = run_mineos(model1,trudata.(itp{1}).periods,'L',ID,0,0,par.inv.verbose);
            end
            if any(strcmp(allpdytp(:,2),'HV')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'HV'),1,'first'));
                [SW.HV.HVr,HVK_new] = run_HVkernel(model1,trudata.(itp{1}).periods,['HV_',ID],1,0,par.inv.verbose);
            end
            for id = 1:length(par.inv.datatypes)
                
                dtype = par.inv.datatypes{id}; pdtyp=parse_dtype(dtype); 
                if ~strcmp(pdtyp{1},'SW'),continue; end
                switch pdtyp{2}
                    case {'Ray','Lov'}
                        swk = predata.(dtype).phV; % record existing phV from kernels
                        predata.(dtype).phV = SW.(pdtyp{2}).(pdtyp{3});
                        swd = predata.(dtype).phV; % record new, precise phV from mineos
                    case 'HV'
                        swk = predata.(dtype).HVr; % record existing phV from kernels
                        predata.(dtype).HVr = SW.(pdtyp{2}).(pdtyp{3});
                        swd = predata.(dtype).HVr; % record new, precise phV from Tanimoto script
                end
                
                if par.inv.verbose
                    fprintf('   %s RMS diff is %.4f\n',par.inv.datatypes{id},rms(swk-swd)); % RMS difference
                end
            end
			newK = true;
		else 
			Ktry = [];
		end
    end % only redo data if model has changed 

%      plot_TRUvsPRE_old(trudata,predata)]

    % continue if any Sp or PS inhomogeneous or nan or weird output
    if ifforwardfail(predata,par)
        fail_chain=fail_chain+1; ifpass=0; 
        fprintf('Forward model error, failchain %.0f\n',fail_chain);  break;
    else
        fail_chain = 0; 
    end

%% =========================  CALCULATE MISFIT  ===========================
    
    % SW weights, if applicable 
    [ SWwt ] = make_SW_weight( par,Kbase,trudata );
    
    [ misfit1 ] = b4_CALC_MISFIT( trudata,predata,par,0,SWwt ); % misfit has structures of summed errors

%% =======================  CALCULATE LIKELIHOOD  =========================
    [ log_likelihood1,misfit1 ] = b5_CALC_LIKELIHOOD( misfit1,trudata,model1.datahparm,par);
    
%     fprintf('MISFITS: Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.SpRF,misfit.PsRF,misfit.SW)
%     fprintf('CHI2S:   Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.chi2_sp,misfit.chi2_ps,misfit.chi2_SW)
    
    fail_chain = 0;
    predat_save1 = predata0;

    end % while ifpass
    
%% ========================  ACCEPTANCE CRITERION  ========================
    [ ifaccept ] = b6_IFACCEPT( log_likelihood1,log_likelihood,temp,p_bd*ifpass);
    
    % ======== PLOT ========  if accept
    if ifaccept && par.inv.verbose && fail_chain==0
        plot_TRUvsPRE( trudata,predata);
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
        
    %% UPDATE KERNEL if needed 
        if newK==true
            Kbase.modelk = model;
            for id = 1:length(par.inv.datatypes)
                dtype = par.inv.datatypes{id};pdtyp = parse_dtype(dtype); 
                if ~strcmp(pdtyp{1},'SW'), continue; end
                if strcmp(pdtyp{2},'HV')
                    K =  HVK_new;   
                else
                    K = run_kernels(trudata.(dtype).periods,pdtyp{2},pdtyp{3},ID,1,0,par.inv.verbose);
                end
                Kbase = populate_Kbase( Kbase,dtype,predata.(dtype).(pdtyp{3}),[],{K} );
            end
            nchain = 0;
        end
    
    else
        if par.inv.verbose, fprintf('  --FAIL--\n'); end
        if newK, delete_mineos_files(ID,'R'); end
        if newK, delete_mineos_files(ID,'L'); end
    end
    
    % restart-chain if immediate failure
    if isinf(log_likelihood), fail_chain=100; break; end 
    
    
%% =========  redo kernel at end of burn in or if chain is too long =======
    if newK == false
        if nchain > par.inv.maxnkchain
            fprintf('\n RECALCULATING KERNEL - too long chain\n'), newK = true;
        elseif ii == par.inv.burnin
            fprintf('\n RECALCULATING KERNEL - end of burn in\n'), newK = true;
        end
        if newK == true
            [Kbase] = make_allkernels(model,Kbase,trudata,ID,par);
            nchain = 0;
        end
    end
    
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

model0_perchain{iii} = model0;
misfits_perchain{iii} = misfits;
allmodels_perchain{iii} = allmodels;
if isfield(trudata,'SW_Ray')
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

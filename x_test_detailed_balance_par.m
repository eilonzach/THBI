clear all
close all


notes = [...
    'Detailed balance test \n',...
    'including p_bd in likelihood term, so possiblity of perturbation fails  \n'...
    'The p_bd is now calculated as just kold/knew, ignoring the q(v1,v0) term.  \n'...
    'Select new knot velocity from gaussian around current value.  \n'...
    'Defining k (N of knots) separately in each layer, rather than overall!  \n'...
        ];
    
% perturb_opt = 'PAPER';
% perturb_opt = 'constant';
% perturb_opt = 'uniformV';
perturb_opt = 'BodinAcc';



%% ------------------------- START ------------------------- 
global projdir THBIpath TRUEmodel
THBIpath = '/Users/zeilon/Documents/MATLAB/BayesianJointInv';
projdir = [THBIpath,'/test_detailed_balance/'];

run([THBIpath,'/a0_STARTUP_BAYES']);
cd(projdir);

%% PARMS
run('bayes_inv_parms')
if strcmp(perturb_opt,'constant')
    par.mod.mantle.kmax = 10;
    par.mod.mantle.kmin = 10;
    par.mod.crust.kmax = 4;
    par.mod.crust.kmin = 4;
end

% par.inv.verbose=false;
par.sta = 'detailed_balance'; par.nwk=[];

STAMP=['DetBa_',perturb_opt]; par.STAMP = STAMP;
resdir = projdir;
mkdir(resdir);
fid = fopen([resdir,'/notes.txt'],'w'); fprintf(fid,notes); fclose(fid);

par_ORIG = par;
% par came from above script
save([resdir,'/par_',STAMP],'par');
% copyfile('bayes_inv_parms.m',[resdir,STAMP,'_parms.m']);

zatdep = [5:5:par.mod.maxz]';

% PRIOR
% parfor kkk = 1:par.inv.nchains
% chainstr = mkchainstr(kkk);
% fprintf('  > Building prior distribution %s from %.0f runs\n',chainstr,min([par.inv.niter,5e5]))
% priors{kkk,1} = a2_BUILD_PRIOR(par,min([par.inv.niter,1e4]),zatdep);
% end
% prior = priors{1};
% fns = fieldnames(prior);
% for kkk = 2:par.inv.nchains 
%     for jj = 1:length(fns)
%         if isscalar(prior.(fns{jj}))
%             prior.(fns{jj}) = prior.(fns{jj}) + priors{kkk}.(fns{jj});
%         elseif strcmp(fns{jj},'zatdep')
%             continue;
%         else
%             prior.(fns{jj}) = [prior.(fns{jj}) ; priors{kkk}.(fns{jj})];
%         end
%     end
% end
% plot_MODEL_SUMMARY(prior,1,[resdir,'/prior_fig.pdf']);
% save([resdir,'/prior'],'prior');
    
%% ---------------------------- INITIATE ----------------------------
profile clear
profile on

% figure(22), clf, set(gcf, 'pos',[480 6 652 799]); hold on, 
% plot(prior.VSmantle(1:150:end,:)',prior.zmantle,'r')
% hp=plot(prior.VSmantle(1:150:end,:)',prior.zmantle,'r');
% set(gca,'ydir','reverse')


%% START DIFFERENT MARKOV CHAINS IN PARALLEL
%% ========================================================================
%% ========================================================================
fprintf('\n =========== STARTING PARALLEL CHAINS ===========\n')
%% ========================================================================
%% ========================================================================
t = now;
parfor iii = 1:par.inv.nchains 
% end
chainstr = mkchainstr(iii);

%% Prep posterior structure
[ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par);

%% Initiate model
ifpass = 0;
% only let starting model out of the loop if it passes conditions
while ifpass==0 % first make sure the starting model satisfies conditions
    model0 = b1_INITIATE_MODEL(par);
    model = model0;
    ifpass = a1_TEST_CONDITIONS( model, par );
end % now we have a starting model!

Nsave = floor(par.inv.niter/par.inv.saveperN);
o = nan(Nsave,1);
% zatdep = linspace(par.mod.sed.hmax+par.mod.crust.hmin+0.1,par.mod.maxz,50)';
dtypes = fieldnames(model(1).datahparm);
posterior = struct('Niter',par.inv.niter*par.inv.nchains,'Nstored',0,...
               'zsed',o,'zmoh',o,...
               'kcrust',o,'kmantle',o,...
               'VSsedtop',o,'VSsedbot',o,'VScrusttop',o,'VScrustbot',o,'VSmanttop',o,...
               'VSmantle',nan(Nsave,length(zatdep)),'zatdep',zatdep,...
               'Zkn_crust',nan(Nsave,par.mod.crust.kmax),'Zkn_mantle',nan(Nsave,par.mod.mantle.kmax),...
               'fdVSmoh',o,'vpvs',o,'cxi',o,'mxi',o,...
               'datahparm',nan(Nsave,length(dtypes)));


%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
ptb = cell({});
acc = zeros(par.inv.niter,1);
crv1 = zeros(par.inv.niter,1);
nchain = 0;
ifaccept=false; 
% reset_likelihood;
misfits.lastlogL = -Inf; 
predata=[];
% not parfor
fprintf('\n =========== STARTING ITERATIONS %s ===========\n',chainstr)
for ii = 1:par.inv.niter
    if rem(ii,100)==0, fprintf('Iteration %s%.0f\n',chainstr,ii); end
    temp = 1;
    
%% ===========================  PERTURB  ===========================  
    if ii==1
        model1 = model; % don't perturb on first run
        p_bd = 1;
        ptb{ii,1} = 'start';
    else
        switch perturb_opt
            case 'constant'
                [model1,ptb{ii,1},p_bd] = b2_PERTURB_MODEL(model,par,temp);
            case 'PAPER'
                [model1,ptb{ii,1},p_bd] = b2_PERTURB_MODEL(model,par,temp);
            case 'BodinAcc'
                [model1,ptb{ii,1},p_bd]  = b2iii_PERTURB_MODEL_BodinAccept(model,par,temp);
            case 'uniformV'
                [model1,ptb{ii,1},p_bd] = b2ii_PERTURB_MODEL_constVbirthp(model,par,temp);
        end  
        
    
%         if strcmp(ptb{ii,1},'Moho_h') && p_bd == 0, [model.crustmparm.knots(end),model1.crustmparm.knots(end)], 
%         end
%         if strcmp(ptb{ii,1},'Moho_h') 
%             junk(size(junk,1)+1,1:3) = [model.crustmparm.knots(end),model1.crustmparm.knots(end),p_bd];
%         end
    
    
        ifpass = a1_TEST_CONDITIONS( model1, par, par.inv.verbose );
        if ifpass==0 && model1.mantmparm.knots(end-1)<par.mod.maxkz
        end
%         ifpass = 1;
        
%         if ~ifpass, 
%         end
    end
    
%         if strcmp(ptb{ii,1},'crust_VS_1') && p_bd == 0, [model.VScrusttop model1.VScrusttop], 
%         end

%     if p_bd ~=0, p_bd = 1; end


%% =====================    ===  ACCEPTANCE CRITERION  ========================
    [ ifaccept ] = b6_IFACCEPT( 0,misfits.lastlogL,temp,p_bd*ifpass );
%     [model.mantmparm.Nkn,model1.mantmparm.Nkn,round(model1.mantmparm.Nkn-model.mantmparm.Nkn),p_bd,ifaccept]
%     if strcmp(ptb{ii,1}(1:end-2),'crust_VS') && ifaccept == 0, [model.crustmparm.VS_sp model1.crustmparm.VS_sp], 
%     end

        
%% ========================  IF ACCEPT ==> STORE  =========================
    if ifaccept 
        model = model1;
        misfits.lastlogL=0;
        acc(ii)=1;
        crv1(ii) = model.crustmparm.VS_sp(1);
    else
        if par.inv.verbose, fprintf('    > fail\n'); end
        acc(ii) = 0;
        crv1(ii) = model.crustmparm.VS_sp(1);
    end
    
    
    % SAVE into 1-D posterior struct
    if mod(ii,par.inv.saveperN)==0
    istor = posterior.Nstored+1;
    posterior.Nstored=istor;
    posterior.ifaccept(istor,1) = ifaccept;
    posterior.zmoh(istor,1) = model.zmoh;
    posterior.fdVSmoh(istor,1) = model.fdVSmoh;
    posterior.kcrust(istor,1) = model.crustmparm.Nkn;
    posterior.kmantle(istor,1) = model.mantmparm.Nkn;
    posterior.VScrusttop(istor,1) = model.crustmparm.VS_sp(1);
    posterior.VScrustbot(istor,1) = model.crustmparm.VS_sp(end);
    posterior.VSmanttop(istor,1) = model.mantmparm.VS_sp(1);
    posterior.VSmantle(istor,:) = linterp(model.z,model.VS,posterior.zatdep);
    posterior.Zkn_crust(istor,1:model.crustmparm.Nkn-2) = model.crustmparm.knots(2:end-1);
    posterior.Zkn_mantle(istor,1:model.mantmparm.Nkn-2) = model.mantmparm.knots(2:end-1);
    posterior.vpvs(istor,1) = model.crustmparm.vpvs;
    posterior.cxi(istor,1) = model.crustmparm.xi;
    posterior.mxi(istor,1) = model.mantmparm.xi;
    for id = 1:length(dtypes)
        posterior.datahparm(istor,id) = model.datahparm.(dtypes{id});
    end
    
    
    end
    
%     if rem(ii,1000)==0 
%         figure(1); clf
% %         for jj = 1:5
%         subplot(1,2,1)
%         histogram(posterior.VScrusttop,[2.475:0.05:4.325]); xlim([2.5 4.3]); hold on
%         subplot(1,2,2)
%         histogram(posterior.VScrustbot,[2.475:0.05:4.325]); xlim([2.5 4.3]); hold on
%             
% %         histogram(posterior.kmantle,[2.5:15.5]); hold on; plot([3:15],posterior.Nstored*(1./[3:15])./sum(1./[3:15]),'r','linewidth',3)
% %         subplot(512);histogram(posterior.kcrust,[1.5:6.5]); hold on; plot([2:6],posterior.Nstored*(1./[2:6])./sum(1./[2:6]),'r','linewidth',3)
%         pause(0.01);
% %     figure(22), hold on
% %     if exist('hp') && iscell(hp)
% %         delete(hp{1}),delete(hp{2}),delete(hp{3}),delete(hp{4}),delete(hp{5}),delete(hp{6}),delete(hp{7}),delete(hp{8})
% %     end
% %     clear hp; hp=cell({});
% %     hp{8}= plot(posterior.VSmantle(1:istor,:)',posterior.zatdep,'--c','linewidth',0.2); 
% %     hp{2} = plot(posterior.VSmantle(max([1,istor-20]):max([1,istor-1]),:)',posterior.zatdep,'b','linewidth',1.8); 
% %     hp{3} = plot(posterior.VSmantle(max([1,istor-40]):max([1,istor-21]),:)',posterior.zatdep,'b','linewidth',1.5); 
% %     hp{4} = plot(posterior.VSmantle(max([1,istor-60]):max([1,istor-41]),:)',posterior.zatdep,'b','linewidth',1.1); 
% %     hp{5} = plot(posterior.VSmantle(max([1,istor-80]):max([1,istor-61]),:)',posterior.zatdep,'b','linewidth',0.8); 
% %     hp{6} = plot(posterior.VSmantle(max([1,istor-100]):max([1,istor-81]),:)',posterior.zatdep,'b','linewidth',0.5); 
% %     hp{7} = plot(posterior.VSmantle(max([1,istor-150]):max([1,istor-101]),:)',posterior.zatdep,'b','linewidth',0.3); 
% %     hp{1}= plot(posterior.VSmantle(istor,:)',posterior.zatdep,'k','linewidth',2); 
%     pause(0.02);
%     end
    
end % on iterations
%% -------------------------- End iteration  ------------------------------
fprintf('\n =========== ENDING ITERATIONS %s ===========\n',char(64+iii))
posteriors{iii} = posterior;
ptbs{iii} = ptb;

end % parfor loop
%% ========================================================================
%% ========================================================================
%% ----------------------- End loop on chains  ----------------------------
%% ========================================================================
%% ========================================================================
fprintf('Duration of entire run: %.0f s\n',(now-t)*86400)
load('prior')

posterior = posteriors{1};
ptb = ptbs{1};
fns = fieldnames(posterior);
for iii = 2:par.inv.nchains 
    for jj = 1:length(fns)
        if isscalar(posterior.(fns{jj}))
            posterior.(fns{jj}) = posterior.(fns{jj}) + posteriors{iii}.(fns{jj});
        elseif strcmp(fns{jj},'zatdep')
            continue;
        else
            posterior.(fns{jj}) = [posterior.(fns{jj}) ; posteriors{iii}.(fns{jj})];
        end
    end
    ptb = {ptb;ptbs{iii}};
end
    
    
plot_db_output
save2pdf(61,['detailed_balance_',perturb_opt]);

return

%% Plot results

figure(91), clf, set(gcf,'pos',[[10 5 1423 800]]);
cls = get(groot,'defaultAxesColorOrder');

%% row 1

axes('pos',[0.05 0.69 0.13 0.24]), cla, hold on
title('Moho depth (km)','fontsize',16)
X = midpts(linspace(par.mod.sed.hmin+par.mod.crust.hmin,par.mod.sed.hmax+par.mod.crust.hmax,20));
No = hist(posterior.zmoh,X)/posterior.Nstored;
Ni = hist(prior.zmoh,X)/prior.Npass;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.sed.hmin+par.mod.crust.hmin par.mod.sed.hmax+par.mod.crust.hmax])

if par.mod.crust.kmin~=par.mod.crust.kmax
axes('pos',[0.21 0.69 0.13 0.24]), cla, hold on
X =par.mod.crust.kmin:par.mod.crust.kmax;
No = hist(posterior.kcrust,X)/posterior.Nstored;
Ni = hist(prior.kcrust,X)/prior.Npass;
predN =  (1./X)./sum(1./X);
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
plot(X,predN,'-b','linewidth',3)
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14), title('Nk - crust','fontsize',16)
end

if par.mod.mantle.kmin~=par.mod.mantle.kmax
axes('pos',[0.37 0.69 0.13 0.24]), cla, hold on
X = par.mod.mantle.kmin:par.mod.mantle.kmax;
No = hist(posterior.kmantle,X)/posterior.Nstored;
Ni = hist(prior.kmantle,X)/prior.Npass;
predN =  (1./X)./sum(1./X);
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
plot(X,predN,'-b','linewidth',3)
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14), title('Nk - mantle','fontsize',16)
end

axes('pos',[0.53 0.69 0.26 0.24]), cla, hold on
X = 0:5:par.mod.maxz;
No = hist([posterior.Zkn_crust(:);posterior.Zkn_mantle(:)],X)/posterior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
No2 = hist(posterior.zmoh(:),X)/posterior.Nstored;
bar(X,No2','facecolor','none','edgecolor',[0.2 0.2 1],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[0 250]), title('knot locations','fontsize',16)

%% row 2

axes('pos',[0.05 0.37 0.13 0.24]), cla, hold on
title('Crust Xi value','fontsize',16)
X = midpts(linspace(par.mod.crust.ximin,par.mod.crust.ximax,20));
Nco = hist(posterior.cxi,X)/posterior.Nstored;
Nci = hist(prior.cxi,X)/prior.Npass;
bar(X,Nco','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Nci','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.crust.ximin par.mod.crust.ximax])

axes('pos',[0.21 0.37 0.13 0.24]), cla, hold on
title('Mantle Xi value','fontsize',16)
X = midpts(linspace(par.mod.mantle.ximin,par.mod.mantle.ximax,20));
Nmo = hist(posterior.mxi,X)/posterior.Nstored;
Nmi = hist(prior.mxi,X)/prior.Npass;
bar(X,Nmo','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Nmi','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.mantle.ximin,par.mod.mantle.ximax])

axes('pos',[0.37 0.37 0.13 0.24]), cla, hold on
X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
No = hist(posterior.VScrusttop,X)/posterior.Nstored;
Ni = hist(prior.VScrusttop,X)/prior.Npass;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.crust.vsmin,par.mod.crust.vsmax]), title('Vs crust top (km/s)','fontsize',16)

axes('pos',[0.53 0.37 0.13 0.24]), cla, hold on
title('Vs crust bot (km/s)','fontsize',16)
X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
No = hist(posterior.VScrustbot,X)/posterior.Nstored;
Ni = hist(prior.VScrustbot,X)/prior.Npass;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.crust.vsmin,par.mod.crust.vsmax])

axes('pos',[0.69 0.37 0.13 0.24]), cla, hold on
X = midpts(linspace(0,30,20));
No = hist(posterior.fdVSmoh,X)/posterior.Nstored;
Ni = hist(prior.fdVSmoh,X)/prior.Npass;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14), title('fractional dVs at moho (%)','fontsize',16)

%% row 3
axes('pos',[0.05 0.06 0.13 0.24]), cla, hold on
title('crust Vp/Vs ratio','fontsize',16)
X = midpts(linspace(par.mod.crust.vpvsmin,par.mod.crust.vpvsmax,20));
No = hist(posterior.vpvs,X)/posterior.Nstored;
Ni = hist(prior.vpvs,X)/prior.Npass;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.crust.vpvsmin,par.mod.crust.vpvsmax]) 


%% velocity at depth
X = midpts(linspace(par.mod.mantle.vsmin,par.mod.mantle.vsmax,20));
zdo = [80,100,120,150,200,250]; for iz = 1:length(zdo), izdo(iz) = crossing(posterior.zatdep,[],zdo(iz));end
for iz = 1:6
axes('pos',[0.10+iz*0.11 0.06 0.08 0.24]), cla, hold on
No = hist(posterior.VSmantle(:,izdo(iz)),X)/posterior.Nstored;
Ni = hist(prior.VSmantle(:,izdo(iz)),X)/prior.Npass;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% legend(num2str(model_summary.zmantle(:)),'location','northwest')
set(gca,'fontsize',14,'xlim',[par.mod.mantle.vsmin,par.mod.mantle.vsmax],'ylim',[0 0.07*diff([par.mod.mantle.vsmin,par.mod.mantle.vsmax])]), 
title(sprintf('Vs at %.0f km',prior.zatdep(izdo(iz))),'fontsize',16)
end

%% hyperparms
for id = 1:length(par.inv.datatypes)
axes('pos',[0.86 (0.06 + (id-1)*0.235) 0.12 0.18]), cla, hold on
X = midpts(linspace(log10(par.mod.data.min_sigma.BW.Ps.def),1,20));
No = hist(log10(posterior.datahparm(:,id)),X)/posterior.Nstored;
% Ni = hist(log10(prior.VSmantle(:,izdo(iz)),X)/prior.Npass;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% % legend(num2str(model_summary.zmantle(:)),'location','northwest')
set(gca,'fontsize',14,'xlim',[-3 0]), title(regexprep(sprintf('%s',par.inv.datatypes{id}),'_','-'),'fontsize',16)
end


%% title
% htit = title_custom([par.sta,' ',par.nwk],0.95,'fontweight','bold','fontsize',25);


%% Save some things
fprintf('  > Saving posterior\n')
save([resdir,'/posterior'],'posterior');
save2pdf(91,['prior2posterior_',perturb_opt]);


return

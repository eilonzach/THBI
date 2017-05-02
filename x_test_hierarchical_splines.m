clear
project_name = 'SYNTHETICS';

run([project_name,'/parms/bayes_inv_parms'])
par.mod.dz=1;
% z0_SYNTH_MODEL_simplemod(par,1);  %close(95)
z0_SYNTH_MODEL_fig2_splines(par,1)
global TRUEmodel
global TLM

%% initiate
par.mod.mantle.kmin = 10;par.mod.mantle.kmax = 10;
model0 = b1_INITIATE_MODEL(par);
model=model0;

par.mod.mantle.kmin = 4;par.mod.mantle.kmax = 20;
par.mod.crust.kmin = 2;par.mod.crust.kmax = 5;

theta = 0.04; Vsmf = inf; Vpmf=inf;

%% loop on iterations
log_like = -inf;
for ii = 1:10000
    fprintf('Iteration %.0f\n',ii); 
    %% perturb
    ifpass = 0;
    % only let perturbed model out of the loop if it passes conditions
    while ifpass==0
        [model1,ptbtype,p_bd ] = b2_PERTURB_MODEL(model,par,1);
        ifpass = a1_TEST_CONDITIONS( model1, par );
        if ~isempty(regexp(ptbtype,'sigma','once')), ifpass = false; end
    end

%     if p_bd~=1, pause,end

    %% misfit
    zzz = 0:0.1:par.mod.maxz; 
    Vsmf1 = norm(linterp(TRUEmodel.z,TRUEmodel.vs,zzz)- linterp(model1.z,model1.VS,zzz));
    Vpmf1 = norm(linterp(TRUEmodel.z,TRUEmodel.vp,zzz)- linterp(model1.z,model1.VP,zzz));
%     Vsmf1 = norm(TRUEmodel.vs - linterp(model1.z,model1.VS,TRUEmodel.z));
    
    %% new likelihood
    log_like1 = -0.5*(Vsmf1.^2 + 0.25*Vpmf1.^2)/theta/theta;
    log_likeR = log_like1 - log_like + log(p_bd)/theta ;
    
    ifaccept=false;
    if log_likeR>0
        ifaccept = true;
    elseif exp(log_likeR) >= random('unif',0,1,1);
        ifaccept = true;
    end
    if ifaccept
        fprintf('> Misfit from %.3f to %.3f\n',Vsmf+Vpmf,Vsmf1+Vpmf1)
        model = model1;
        log_like = log_like1;
        Vsmf  = Vsmf1;
        Vpmf  = Vpmf1;
        plot_MOD_TRUEvsTRIAL(TRUEmodel,model);
        pause(0.01)
    end
end
%% see parameterisation
plot_PARAMETERISATION(model)

%% see if layerised fit good too
[zlayt,zlayb,Vslay,Vplay,rholay] = ...
    layerise(model.z,model.VS,par.forc.mindV,0,model.VP,model.rho); 

figure(44), clf, set(gcf,'pos',[680 436 705 662])
subplot(121), hold on
plot(TRUEmodel.VS,TRUEmodel.z,'--g','linewidth',2)
plot(reshape([TLM.Vs';TLM.Vs'],2*TLM.nlay,1),reshape([TLM.zlayt';TLM.zlayb'],2*TLM.nlay,1),'-ok')
plot(reshape([Vslay';Vslay'],2*length(zlayt),1),reshape([zlayt';zlayb'],2*length(zlayt),1),'-or','linewidth',1.5)
set(gca,'ydir','reverse','fontsize',18); xlabel('VS','fontsize',20),ylabel('Depth (km)','fontsize',20)
subplot(122), hold on
plot(TRUEmodel.VP,TRUEmodel.z,'--g','linewidth',2)
plot(reshape([TLM.Vp';TLM.Vp'],2*TLM.nlay,1),reshape([TLM.zlayt';TLM.zlayb'],2*TLM.nlay,1),'-ok')
plot(reshape([Vplay';Vplay'],2*length(zlayt),1),reshape([zlayt';zlayb'],2*length(zlayt),1),'-or','linewidth',1.5)
set(gca,'ydir','reverse','fontsize',18); xlabel('VS','fontsize',20)

%% see what data look like
% true
[ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
% model
phV = run_mineos(model,trudata.SW.periods,'start',0,0,1);
K   = run_kernels(model,trudata.SW.periods,'start',1,0,1);
Kbase = struct('modelk',model,'phV',phV,'K',{K});
predata = b3_FORWARD_MODEL( model,Kbase,par,trudata,'eg',0);
predata = rmfield(predata,{'PsRF_lo','SpRF_lo'});
plot_TRUvsPRE(trudata,predata)

%% see how similar final splines are
figure(46), clf, set(gcf,'pos',[100 100 400 800])
ax = gca;
cla(ax);
hold(ax,'on')
% crust splines
for is = 1:size(model.crustmparm.splines,2)
    plot(ax,model.crustmparm.splines(:,is)*model.crustmparm.VS_sp(is),model.crustmparm.z_sp,'-g','linewidth',3)
end
for is = 1:size(TRUEmodel.crustmparm.splines,2)
    plot(ax,TRUEmodel.crustmparm.splines(:,is)*TRUEmodel.crustmparm.VS_sp(is),TRUEmodel.crustmparm.z_sp,':k','linewidth',1.5)
end
% mantle splines
for is = 1:size(model.mantmparm.splines,2)
    plot(ax,model.mantmparm.splines(:,is)*model.mantmparm.VS_sp(is),model.mantmparm.z_sp,'-b','linewidth',3)
end
for is = 1:size(TRUEmodel.mantmparm.splines,2)
    plot(ax,TRUEmodel.mantmparm.splines(:,is)*TRUEmodel.mantmparm.VS_sp(is),TRUEmodel.mantmparm.z_sp,':k','linewidth',1.5)
end

% BOUNDARIES
% crust/mantle boundary
plot(ax,[0,ax.XLim(2)],max(model.zmoh)*[1,1],'-','linewidth',3,'color',[0.5 0.5 0.5])
plot(ax,[0,ax.XLim(2)],max(TRUEmodel.zmoh)*[1,1],':k','linewidth',1)
text(ax,0.5,max(model.zmoh),'Moho','HorizontalAlignment','center','verticalalignment','top','fontsize',15,'fontweight','bold','color',[0.5 0.5 0.5])


set(ax,'ydir','reverse','fontsize',14,'box','on','layer','top','linewidth',3,'xtick',[0 1])
title(ax,'Model parameterisation','fontsize',20)
xlabel(ax,'Weight','fontsize',17,'fontweight','bold')
ylabel(ax,'Depth (km)','fontsize',17,'fontweight','bold')



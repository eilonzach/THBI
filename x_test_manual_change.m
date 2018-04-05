%% take current model, perturb it manually, see if fits improve
% profile clear
% profile on
% close

%% current model data and misfit
% [ model ] = make_mod_from_parms( model,par );
cdata = b3_FORWARD_MODEL( model,Kbase,par,trudata,id,0);
% window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ trudata ] = predat_process( trudata,dtype,par);
end
% [cdata.SW.phV] = run_mineos(model,trudata.SW.periods,'cmodmake');
[ cmisfit ] = b4_CALC_MISFIT( trudata,cdata,par,0,SWwt ); % misfit has structures of summed errors
[ clog_likelihood,cmisfit ] = b5_CALC_LIKELIHOOD( cmisfit,trudata,model.datahparm,par );
    
%% perturb model
model2 = model;

% model2 = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','crust','remove'}, 3); 
% [model2] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',1}, -0.2);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','crust','vpvs',4}, -0.0183);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','add'},[90 +0.2]);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','add'},[95 +0.1]);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','add'},[98 +0.04]);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','add'},[110 -0.2]);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11); 
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11); 
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11);
[model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',8},-0.1);
[model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',9},-0.1);
[model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',10},-0.1);
[model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',12},+0.2);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','crust','V',2},-0.05);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','crust','V',3},-0.15);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',1},-0.4);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},8); 
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',5},-0.1);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',6},-0.05);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',7},-0.2);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',9},+0.1);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',10},+0.17);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},3);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','zkn',3},-13);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',7},-0.05);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',10},+0.05);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',2},-0.01);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','sed',1},+0.12);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','crust','vpvs'},+0.11);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_disc','moh','dV'},-0.5);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_disc','moh','h'},-4);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_disc','sed','h'},-0.1);
% model2 = TRUEmodel;
% model2.datahparm = model.datahparm;
% model2.crustmparm = TRUEmodel.crustmparm;

[ model2 ] = make_mod_from_parms( model2,par );


% compute perturbation
modptb = calc_Vperturbation(model,model2);
ptbnorm = norm(modptb.dvsv) + norm(modptb.dvsh);

%% perturbed model data and misfit
predata = b3_FORWARD_MODEL( model2,Kbase,par,trudata,id,0);
% window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ trudata ] = predat_process( trudata,dtype,par);
end
% [predata.SW.phV] = run_mineos(model2,trudata.SW.periods,'cmodmake');
[ misfit ] = b4_CALC_MISFIT( trudata,predata,par,0,SWwt ); % misfit has structures of summed errors
[ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,trudata,model2.datahparm,par );


% plot true, current, and trial models
plot_MOD_TRUEvsTRIAL( TRUEmodel, model )
subplot(1,2,1), hold on, h(1) = plot(model2.VS,model2.z,'b--');
subplot(1,2,2), hold on, h(2) = plot(model2.VP,model2.z,'b--');

%% plot true, current, and trial data
axs = plot_TRUvsPRE_WAVEFORMS( trudata,cdata); 

% axes(axs(1)), hold on,plot(predata.BW_Ps_lo.tt,predata.BW_Ps_lo.PSV(:,1)./maxab(predata.BW_Ps_lo.PSV(:)),'c','linewidth',2),
% axes(axs(2)), hold on,plot(predata.BW_Ps_lo.tt,predata.BW_Ps_lo.PSV(:,2)./maxab(predata.BW_Ps_lo.PSV(:)),'c','linewidth',2),
% axes(axs(3)), hold on,plot(predata.BW_Sp_lo.tt,predata.BW_Sp_lo.PSV(:,1)./maxab(predata.BW_Sp_lo.PSV(:)),'c','linewidth',2),
% axes(axs(4)), hold on,plot(predata.BW_Sp_lo.tt,predata.BW_Sp_lo.PSV(:,2)./maxab(predata.BW_Sp_lo.PSV(:)),'c','linewidth',2),
for ipp = 1:length(predata.BW_Ps)
    axes(axs(5)), hold on,plot(predata.BW_Ps(ipp).tt,predata.BW_Ps(ipp).PSV(:,1)./maxab(predata.BW_Ps(ipp).PSV(:)),'c','linewidth',2),
    axes(axs(6)), hold on,plot(predata.BW_Ps(ipp).tt,predata.BW_Ps(ipp).PSV(:,2)./maxab(predata.BW_Ps(ipp).PSV(:)),'c','linewidth',2),
end
for iss = 1:length(predata.BW_Sp)
    axes(axs(7)), hold on,plot(predata.BW_Sp(iss).tt,predata.BW_Sp(iss).PSV(:,1)./maxab(predata.BW_Sp(iss).PSV(:)),'c','linewidth',2),
    axes(axs(8)), hold on,plot(predata.BW_Sp(iss).tt,predata.BW_Sp(iss).PSV(:,2)./maxab(predata.BW_Sp(iss).PSV(:)),'c','linewidth',2),
end
axes(axs(9)), hold on,plot(predata.SW_Ray_phV.periods,predata.SW_Ray_phV.phV,'c.-','linewidth',3,'markersize',40);

pause(0.01)
figure(85)

%% report changes in misfit and likelihood
fprintf('=============================\n')
fprintf('              model   model1\n')
fprintf('-----------------------------\n')
fprintf('BW_Psmisfit %8.3f  %8.3f  %8.3f\n',mean(cmisfit.E2.BW_Ps),mean(misfit.E2.BW_Ps),mean(cmisfit.E2.BW_Ps)-mean(misfit.E2.BW_Ps))
fprintf('BW_Spmisfit %8.3f  %8.3f  %8.3f\n',mean(cmisfit.E2.BW_Sp),mean(misfit.E2.BW_Sp),mean(cmisfit.E2.BW_Sp)-mean(misfit.E2.BW_Sp))
% fprintf('BW_Psmisflo %8.3f  %8.3f  %8.3f\n',cmisfit.BW_Ps_lo,misfit.BW_Ps_lo,cmisfit.BW_Ps_lo-misfit.BW_Ps_lo)
% fprintf('BW_Spmisflo %8.3f  %8.3f  %8.3f\n',cmisfit.BW_Sp_lo,misfit.BW_Sp_lo,cmisfit.BW_Sp_lo-misfit.BW_Sp_lo)
fprintf('SWmisfit    %8.3f  %8.3f  %8.3f\n',cmisfit.E2.SW_Ray_phV,misfit.E2.SW_Ray_phV,cmisfit.E2.SW_Ray_phV-misfit.E2.SW_Ray_phV)
fprintf('Ps chi2     %8.3f  %8.3f  %8.3f\n',mean(cmisfit.chi2.BW_Ps),mean(misfit.chi2.BW_Ps),mean(cmisfit.chi2.BW_Ps)-mean(misfit.chi2.BW_Ps))
fprintf('Sp chi2     %8.3f  %8.3f  %8.3f\n',mean(cmisfit.chi2.BW_Sp),mean(misfit.chi2.BW_Sp),mean(cmisfit.chi2.BW_Sp)-mean(misfit.chi2.BW_Sp))
% fprintf('Ps_lo chi2 %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_ps_lo,misfit.chi2_ps_lo,cmisfit.chi2_ps_lo-misfit.chi2_ps_lo)
% fprintf('Sp_lo chi2 %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_sp_lo,misfit.chi2_sp_lo,cmisfit.chi2_sp_lo-misfit.chi2_sp_lo)
fprintf('SW chi2     %8.3f  %8.3f  %8.3f\n',cmisfit.chi2.SW_Ray_phV,misfit.chi2.SW_Ray_phV,cmisfit.chi2.SW_Ray_phV-misfit.chi2.SW_Ray_phV)
fprintf('Like''hood   %8.3f  %8.3f  %8.3f\n',clog_likelihood,log_likelihood,log_likelihood-clog_likelihood)
% 
% profile viewer

% figure(11)
% subplot(1,3,1),ylim([0,50]),xlim([3.4,4.4]);
% subplot(1,3,2),ylim([0,50]),xlim([6,8]);
% subplot(1,3,3),ylim([0,50]),xlim([2.6,3.4])
return
%% loop over hyperparamters
sigs = logspace(-3.5,0,60)';
l_hoods = zeros(length(sigs),3);
model1.datahparm = model.datahparm;
for iips = 1:length(sigs)
model1.datahparm.sigmaBW_Ps = sigs(iips);
[ l_hoods(iips,1) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaBW_Ps = model.datahparm.sigmaBW_Ps;
for iisp = 1:length(sigs)
model1.datahparm.sigmaBW_Sp = sigs(iisp);
[ l_hoods(iisp,2) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaBW_Sp = model.datahparm.sigmaBW_Sp;
for iisw = 1:length(sigs)
model1.datahparm.sigmaSW = sigs(iisw);
[ l_hoods(iisw,3) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaSW = model.datahparm.sigmaSW;
figure(45),clf
loglog(sigs,l_hoods);


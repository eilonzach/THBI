%% take current model, perturb it manually, see if fits improve
% profile clear
% profile on

%% current model data and misfit
% [ model ] = make_mod_from_parms( model,par );
cdata = b3_FORWARD_MODEL( model,Kbase,par,trudata,id,0);
cdata = predat_process(cdata,'PsRF',par);
cdata = predat_process(cdata,'SpRF',par);
cdata = predat_process(cdata,'PsRF_lo',par);
cdata = predat_process(cdata,'SpRF_lo',par);
% [cdata.SW.phV] = run_mineos(model,trudata.SW.periods,'cmodmake');
[ cmisfit ] = b4_CALC_MISFIT( trudata,cdata,par,0,SWwt ); % misfit has structures of summed errors
[ clog_likelihood,cmisfit ] = b5_CALC_LIKELIHOOD( cmisfit,trudata,model.datahparm,par );
    
%% perturb model
model2 = model;

% model2 = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','crust','remove'}, 3); 
[model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','crust','V',2}, -0.2);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','crust','vpvs',4}, -0.0183);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11); 
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11); 
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},11);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',11},+0.5);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},8); 
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',5},-0.1);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',4},-0.1);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Mmod','mantle','remove'},3);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','zkn',3},-13);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',7},-0.05);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',10},+0.05);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',2},-0.01);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','mantle','V',1},-0.02);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_Vmod','crust','V',4},+0.02);
% [model2,ptb] = x_CUSTOM_PERTURB_MODEL( model2,par, {'ptb_disc','moh','dV'},-0.04);
% model2 = TRUEmodel;
% model2.datahparm = model.datahparm;
% model2.crustmparm = TRUEmodel.crustmparm;

[ model2 ] = make_mod_from_parms( model2,par );


% compute perturbation
modptb = calc_Vperturbation(model,model2);
ptbnorm = norm(modptb.dvsv) + norm(modptb.dvsh);

%% perturbed model data and misfit
predata = b3_FORWARD_MODEL( model2,Kbase,par,trudata,id,0);
predata = predat_process(predata,'PsRF',par);
predata = predat_process(predata,'SpRF',par);
predata = predat_process(predata,'PsRF_lo',par);
predata = predat_process(predata,'SpRF_lo',par);
% [predata.SW.phV] = run_mineos(model2,trudata.SW.periods,'cmodmake');
[ misfit ] = b4_CALC_MISFIT( trudata,predata,par,0,SWwt ); % misfit has structures of summed errors
[ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,trudata,model2.datahparm,par );


%% plot true, current, and trial models
plot_MOD_TRUEvsTRIAL( TRUEmodel, model )
subplot(1,2,1), hold on, h(1) = plot(model2.VS,model2.z,'b--');
subplot(1,2,2), hold on, h(2) = plot(model2.VP,model2.z,'b--');

%% plot true, current, and trial data
axs = plot_TRUvsPRE_WAVEFORMS( trudata,cdata); 

axes(axs(1)), hold on,plot(predata.PsRF_lo.tt,predata.PsRF_lo.ZRT(:,1),'c','linewidth',2),
axes(axs(2)), hold on,plot(predata.PsRF_lo.tt,predata.PsRF_lo.ZRT(:,2),'c','linewidth',2),
axes(axs(3)), hold on,plot(predata.SpRF_lo.tt,predata.SpRF_lo.ZRT(:,1),'c','linewidth',2),
axes(axs(4)), hold on,plot(predata.SpRF_lo.tt,predata.SpRF_lo.ZRT(:,2),'c','linewidth',2),
axes(axs(5)), hold on,plot(predata.PsRF.tt,predata.PsRF.ZRT(:,1),'c','linewidth',2),
axes(axs(6)), hold on,plot(predata.PsRF.tt,predata.PsRF.ZRT(:,2),'c','linewidth',2),
axes(axs(7)), hold on,plot(predata.SpRF.tt,predata.SpRF.ZRT(:,1),'c','linewidth',2),
axes(axs(8)), hold on,plot(predata.SpRF.tt,predata.SpRF.ZRT(:,2),'c','linewidth',2),
axes(axs(9)), hold on,plot(predata.SW.periods,predata.SW.phV,'c.-','linewidth',3,'markersize',40);

pause(0.01)
figure(85)

%% report changes in misfit and likelihood
fprintf('=============================\n')
fprintf('              model   model1\n')
fprintf('-----------------------------\n')
fprintf('PsRFmisfit %8.3f  %8.3f   %8.3f\n',cmisfit.PsRF,misfit.PsRF,cmisfit.PsRF-misfit.PsRF)
fprintf('SpRFmisfit %8.3f  %8.3f   %8.3f\n',cmisfit.SpRF,misfit.SpRF,cmisfit.SpRF-misfit.SpRF)
fprintf('PsRFmisflo %8.3f  %8.3f   %8.3f\n',cmisfit.PsRF_lo,misfit.PsRF_lo,cmisfit.PsRF_lo-misfit.PsRF_lo)
fprintf('SpRFmisflo %8.3f  %8.3f   %8.3f\n',cmisfit.SpRF_lo,misfit.SpRF_lo,cmisfit.SpRF_lo-misfit.SpRF_lo)
fprintf('SWmisfit   %8.3f  %8.3f   %8.3f\n',cmisfit.SW,misfit.SW,cmisfit.SW-misfit.SW)
fprintf('Ps chi2    %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_ps,misfit.chi2_ps,cmisfit.chi2_ps-misfit.chi2_ps)
fprintf('Sp chi2    %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_sp,misfit.chi2_sp,cmisfit.chi2_sp-misfit.chi2_sp)
fprintf('Ps_lo chi2 %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_ps_lo,misfit.chi2_ps_lo,cmisfit.chi2_ps_lo-misfit.chi2_ps_lo)
fprintf('Sp_lo chi2 %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_sp_lo,misfit.chi2_sp_lo,cmisfit.chi2_sp_lo-misfit.chi2_sp_lo)
fprintf('SW chi2    %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_SW,misfit.chi2_SW,cmisfit.chi2_SW-misfit.chi2_SW)
fprintf('Like''hood  %8.3f  %8.3f   %8.3f\n',clog_likelihood,log_likelihood,log_likelihood-clog_likelihood)
% 
% profile viewer

figure(11)
subplot(1,3,1),ylim([0,50]),xlim([3.4,4.4]);subplot(1,3,2),ylim([0,50]),xlim([6,8]);subplot(1,3,3),ylim([0,50]),xlim([2.6,3.4])
return
%% loop over hyperparamters
sigs = logspace(-3.5,0,60)';
l_hoods = zeros(length(sigs),3);
model1.datahparm = model.datahparm;
for iips = 1:length(sigs)
model1.datahparm.sigmaPsRF = sigs(iips);
[ l_hoods(iips,1) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaPsRF = model.datahparm.sigmaPsRF;
for iisp = 1:length(sigs)
model1.datahparm.sigmaSpRF = sigs(iisp);
[ l_hoods(iisp,2) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaSpRF = model.datahparm.sigmaSpRF;
for iisw = 1:length(sigs)
model1.datahparm.sigmaSW = sigs(iisw);
[ l_hoods(iisw,3) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaSW = model.datahparm.sigmaSW;
figure(45),clf
loglog(sigs,l_hoods);


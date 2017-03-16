%% take current trial model compare apparent misfit to real misfit

%confirm data
tdata = trudata;
[tdata.SW.phV] = run_mineos(TRUEmodel,trudata.SW.periods,1,'tm',0);

% model1 = model;
% model1.mantmparm.VS_kn(end) = TRUEmodel.mantmparm.VS_kn(end);
% [ model1 ] = make_mod_from_parms( model1,par );
% 

adata = b3_FORWARD_MODEL( model1,Kbase,par,trudata,'kma',0);
% [cdata.SW.phV,~,Ktry] = run_mineos(model,trudata.SW.periods,1,'cmodmake',0);
[ amisfit ] = b4_CALC_MISFIT( tdata,adata,par,0 ); % misfit has structures of summed errors
[ alog_likelihood,amisfit ] = b5_CALC_LIKELIHOOD( amisfit,tdata,model1.datahparm,par );
  


bdata = b3_FORWARD_MODEL( model1,Kbase,par,trudata,'kmb',0);
[bdata.SW.phV] = run_mineos( model1,trudata.SW.periods,1,'kmb',0);
% [cdata.SW.phV,~,Ktry] = run_mineos(model,trudata.SW.periods,1,'cmodmake',0);
[ bmisfit ] = b4_CALC_MISFIT( tdata,bdata,par,0 ); % misfit has structures of summed errors
[ blog_likelihood,bmisfit ] = b5_CALC_LIKELIHOOD( bmisfit,tdata,model1.datahparm,par );
 
%% plots
figure(34), clf, set(gcf,'pos',[380 410 1424 688])
% dispersion curves
subplot(3,4,[1,2,5,6]);hold on, set(gca,'fontsize',14)
h(1) = plot(trudata.SW.periods,trudata.SW.phV,'m.-','linewidth',3,'markersize',40);
h(2) = plot(adata.SW.periods,adata.SW.phV,'b.-','linewidth',3,'markersize',40);
h(3) = plot(bdata.SW.periods,bdata.SW.phV,'c.-','linewidth',3,'markersize',40);
hl = legend(h,{'TRUE','kernels','MINEOS'},'location','northwest'); set(hl,'fontsize',18);
xlabel('Period (s)','fontsize',16)
ylabel('Phase velocity (km/s)','fontsize',16)
% xlim([60 200])
% ylim([3.8 4.5])
subplot(3,4,[9,10]);hold on, set(gca,'fontsize',14), 
h(1) = plot(trudata.SW.periods,trudata.SW.phV,'m.-','linewidth',3,'markersize',40);
h(2) = plot(adata.SW.periods,adata.SW.phV,'b.-','linewidth',3,'markersize',40);
h(3) = plot(bdata.SW.periods,bdata.SW.phV,'c.-','linewidth',3,'markersize',40);
xlim([60 200])
ylim([3.8 4.6])
% VS
subplot(3,4,[3,7,11]), hold on;
plot(TRUEmodel.VS,TRUEmodel.z,'-k','Linewidth',2);
plot(model1.VS,model1.z,'-r','Linewidth',2);
set(gca,'ydir','reverse','fontsize',14,'xlim',[1.5 5]);
title('Vs','fontsize',20)
xlabel('Vs (km/s)','fontsize',16)
ylabel('Depth (km)','fontsize',16)

subplot(3,4,[4,8,12]), hold on;
plot(TRUEmodel.VP,TRUEmodel.z,'-k','Linewidth',2);
plot(model1.VP,model1.z,'-r','Linewidth',2);
set(gca,'ydir','reverse','fontsize',14,'xlim',[3.5 9]);
title('Vp','fontsize',20)
xlabel('Vp (km/s)','fontsize',16)
ylabel('Depth (km)','fontsize',16)
            
fprintf('    RMS diff is %.4f\n',rms(adata.SW.phV-bdata.SW.phV)); % RMS difference



%% report changes in misfit and likelihood
fprintf('=============================\n')
fprintf('             kernel    MINEOS\n')
fprintf('-----------------------------\n')
fprintf('PsRFmisfit %8.3f  %8.3f   %8.3f\n',amisfit.PsRF,bmisfit.PsRF,amisfit.PsRF-misfit.PsRF)
fprintf('SpRFmisfit %8.3f  %8.3f   %8.3f\n',amisfit.SpRF,bmisfit.SpRF,amisfit.SpRF-misfit.SpRF)
fprintf('SWmisfit   %8.3f  %8.3f   %8.3f\n',amisfit.SW,bmisfit.SW,amisfit.SW-misfit.SW)
fprintf('Ps chi2    %8.3f  %8.3f   %8.3f\n',amisfit.chi2_ps,bmisfit.chi2_ps,amisfit.chi2_ps-bmisfit.chi2_ps)
fprintf('Sp chi2    %8.3f  %8.3f   %8.3f\n',amisfit.chi2_sp,bmisfit.chi2_sp,amisfit.chi2_sp-bmisfit.chi2_sp)
fprintf('SW chi2    %8.3f  %8.3f   %8.3f\n',amisfit.chi2_SW,bmisfit.chi2_SW,amisfit.chi2_SW-bmisfit.chi2_SW)
fprintf('Like''hood  %8.3f  %8.3f   %8.3f\n',alog_likelihood,blog_likelihood,blog_likelihood-alog_likelihood)
% 

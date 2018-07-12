%% take current trial model compare apparent misfit to real misfit

%confirm data
% must be in this directory to run
clear 
addpath('/Users/zeilon/Work/codes/HV_Tanimoto/matlab_to_HVkernel');
THBIpath = '/Users/zeilon/Documents/MATLAB/BayesianJointInv';
run([THBIpath,'/a0_STARTUP_BAYES']);
cd(THBIpath);
run test_HVkernel/parms_HVkernel_test.m


swperiods = [16,20,24,28,32,36,40,50,60,70,80,90]'; 
Np = length(swperiods);
niter = 1e2;
nchain = 1;

for jj = 1:nchain
%% Starting model, used for kernels
chainstr = mkchainstr(jj);
ifpass = 0;
while ifpass==0 % first make sure the starting model satisfies conditions
    model0 = b1_INITIATE_MODEL(par);
    ifpass = a1_TEST_CONDITIONS( model0, par );
    try
        [HVr0,HVK0,phV0] = run_HVkernel(model0,swperiods,['test',chainstr],1,0,par.inv.verbose);
    catch
        ifpass = 0;
    end
end
model = model0;

%% perturb model
HVr_true = zeros(niter,Np);
HVr_kern = zeros(niter,Np);
ptbnorm = zeros(niter,1);

for ii = 1:niter
    if rem(ii,50)==0 || ii==1, fprintf('%s%.0f\n',chainstr,ii);end
	temp = 1;
    ifpass = 0;
    itry = 0;
    while ifpass == 0 && itry < 10
        itry = itry+1;
        try
            [model1,~,p_bd] = b2_PERTURB_MODEL(model,par,temp);
            ifpass = a1_TEST_CONDITIONS( model1, par, par.inv.verbose  );
            if ~ifpass, continue; end

            [ modptb ] = calc_Vperturbation( model0,model1);
            HVr_true(ii,:) = run_HVkernel(model1,swperiods,['test',chainstr,num2str(ii)],1,0,par.inv.verbose);
        catch
            ifpass = 0;
        end
    end
    model = model1;
	ptbnorm(ii) = norm(modptb.dvsav) + norm(modptb.dvpav);
    
    dvs_vs = linterp(modptb.Z,modptb.dvsav,0.5*(HVK0{1}.Z1+HVK0{1}.Z2));
    dvp_vp = linterp(modptb.Z,modptb.dvpav,0.5*(HVK0{1}.Z1+HVK0{1}.Z2));
    drh_rh = linterp(modptb.Z,modptb.drho,0.5*(HVK0{1}.Z1+HVK0{1}.Z2));
    
    for ik = 1:length(HVK0)
        HVr_kern(ii,ik) = HVr0(ik) - sum(dvs_vs.*HVK0{ik}.Kzh_Vs + dvp_vp.*HVK0{ik}.Kzh_Vp + drh_rh.*HVK0{ik}.Kzh_rho);
    end
%     figure(6);clf;plot([HVr0,HVr_true(ii,:)',HVr_kern(ii,:)'],'-o')
end


HVr_true_s{jj} = HVr_true;
HVr_kern_s{jj} = HVr_kern;
ptbnorm_s{jj} = ptbnorm;

end

% save('test_HVkernel/results2','HVr_true_s','HVr_kern_s','ptbnorm_s','swperiods')

figure(5); clf, hold on;
ip = randi(Np); 
ip = 7; swperiods(ip)
for jj = 1:nchain-3
HVR_error = 100*(HVr_true_s{jj}(:,ip)-HVr_kern_s{jj}(:,ip))./HVr_true_s{jj}(:,ip);
plot(ptbnorm_s{jj},abs(HVR_error))
end



return

% HVtool

    model = model0;


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

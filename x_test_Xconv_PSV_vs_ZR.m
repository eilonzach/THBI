% script to test the minimising P-SV correlation technique of estimating
% surface velocity
clear 
close all
%% parms
inc = 20;
synthperiod = 3;
samprate = 20;
nsamps = 2^12;
wdo = [-30 30];
P_err_upwt = 1;
noisestd = .05;

%% TRUE model
moh = 40;
Z  = [linspace(0,moh,10) linspace(moh,250,25) 260]';
Vs = [linspace(3.3,3.4,10) linspace(3.9,4.2,25) 4.2]';
Vp = 1.8*Vs;
rho = mantle_vs2rho(Vs,Z);
[zlayt,zlayb,Vs_lay,Vp_lay,rho_lay] = layerise(Z,Vs,0.01,0,Vp,rho);
laymod = struct('nlay',length(Vs_lay),'zlayt',zlayt,'zlayb',zlayb,'Vp',Vp_lay,'Vs',Vs_lay,'rho',rho_lay);
raypp = sind(inc)/Vp(end);
rayps = sind(inc)/Vs(end);

%% 'TRUE' Ps synth
[dat_Ps,tt_p] = run_propmat(laymod,'test1','Ps',samprate,inc,synthperiod,nsamps);
dat_Ps = dat_Ps(:,[3,1,2]); % in Z,R,T
tt01 = mean(tt_p(dat_Ps(:,1)==min(dat_Ps(:,1))));% estimate main P-wave arrival time from first big downswing
tt_p = tt_p - tt01; tt_p = round_level(tt_p,0.001);
inwind = (tt_p >= wdo(1)) & (tt_p < wdo(2)); 
dat_Ps = dat_Ps(inwind,:);
tt_p = tt_p(inwind);
true_dat_Ps = flat_hanning_win(tt_p,dat_Ps,wdo(1),wdo(2),3); % 3s taper
true_dat_Ps = true_dat_Ps/max(max(abs(true_dat_Ps)));
% transform to P-SV
[true_dat_Ps_PSV(:,1),true_dat_Ps_PSV(:,2)] = Rotate_XZ_to_PSV(dat_Ps(:,2),-dat_Ps(:,1),Vp(1),Vs(1),raypp);
true_dat_Ps_PSV = true_dat_Ps_PSV/max(max(abs(true_dat_Ps_PSV)));
% noisify
true_dat_Ps = true_dat_Ps + random('norm',0,noisestd,size(true_dat_Ps));
true_dat_Ps_PSV = true_dat_Ps_PSV + random('norm',0,noisestd,size(true_dat_Ps_PSV));
% norms
tru_nrm_ZR_ps  = norm(true_dat_Ps(:,1)) + norm(true_dat_Ps(:,2));
tru_nrm_PSV_ps = norm(true_dat_Ps_PSV(:,1)) + norm(true_dat_Ps_PSV(:,2));


%% 'TRUE' Sp synth
[dat_Sp,tt_s] = run_propmat(laymod,'test1','Sp',samprate,inc,synthperiod,nsamps);

dat_Sp = dat_Sp(:,[3,1,2]); % in Z,R,T
tt01 = mean(tt_s(dat_Sp(:,1)==min(dat_Sp(:,1))));% estimate main P-wave arrival time from first big downswing
tt_s = tt_s - tt01; tt_s = round_level(tt_s,0.001);
inwind = (tt_s >= wdo(1)) & (tt_s < wdo(2)); 
dat_Sp = dat_Sp(inwind,:);
tt_s = tt_s(inwind);
true_dat_Sp = flat_hanning_win(tt_s,dat_Sp,wdo(1),wdo(2),3); % 3s taper
true_dat_Sp = true_dat_Sp/max(max(abs(true_dat_Sp)));
% transform to P-SV
[true_dat_Sp_PSV(:,1),true_dat_Sp_PSV(:,2)] = Rotate_XZ_to_PSV(dat_Sp(:,2),-dat_Sp(:,1),Vp(1),Vs(1),rayps);
true_dat_Sp_PSV = true_dat_Sp_PSV/max(max(abs(true_dat_Sp_PSV)));
% noisify
true_dat_Sp = true_dat_Sp+ random('norm',0,noisestd,size(true_dat_Sp));
true_dat_Sp_PSV = true_dat_Sp_PSV + random('norm',0,noisestd,size(true_dat_Sp_PSV));
% norms
tru_nrm_ZR_sp  = norm(true_dat_Sp(:,1)) + norm(true_dat_Sp(:,2));
tru_nrm_PSV_sp = norm(true_dat_Sp_PSV(:,1)) + norm(true_dat_Sp_PSV(:,2));

%% Grid search through models
moh_try = [20:5:55]';
vpvs_try = [1.5:0.05:2.0]';
misfit_PSV = zeros(length(moh_try),length(vpvs_try),2);
misfit_ZR = zeros(length(moh_try),length(vpvs_try),2);
for imoh = 1:length(moh_try)
for ivc = 1:length(vpvs_try)
    %% FORWARD MODEL
    % TRIAL model
    Z  = [linspace(0,moh_try(imoh),10) linspace(moh_try(imoh),250,25) 260]';
    Vs = [linspace(3.3,3.4,10) linspace(3.9,4.2,25) 4.2]';
    Vp = vpvs_try(ivc)*Vs;
    rho = mantle_vs2rho(Vs,Z);
    [zlayt,zlayb,Vs_lay,Vp_lay,rho_lay] = layerise(Z,Vs,0.01,0,Vp,rho);
    laymod = struct('nlay',length(Vs_lay),'zlayt',zlayt,'zlayb',zlayb,'Vp',Vp_lay,'Vs',Vs_lay,'rho',rho_lay);
    % TRIAL Ps synth
    [dat_Ps,tt_p] = run_propmat(laymod,'test1','Ps',samprate,inc,synthperiod,nsamps);
    dat_Ps = dat_Ps(:,[3,1,2]); % in Z,R,T
    tt01 = mean(tt_p(dat_Ps(:,1)==min(dat_Ps(:,1))));% estimate main P-wave arrival time from first big downswing
    tt_p = tt_p - tt01; tt_p = round_level(tt_p,0.001);
    inwind = (tt_p >= wdo(1)) & (tt_p < wdo(2)); 
    dat_Ps = dat_Ps(inwind,:);
    tt_p = tt_p(inwind);
    try_dat_Ps = flat_hanning_win(tt_p,dat_Ps,wdo(1),wdo(2),3); % 3s taper
    try_dat_Ps = try_dat_Ps/max(max(abs(try_dat_Ps)));
    try_nrm_ZR_ps = norm(try_dat_Ps(:,1)) + norm(try_dat_Ps(:,2));
    % TRIAL Sp synth
    [dat_Sp,tt_s] = run_propmat(laymod,'test1','Sp',samprate,inc,synthperiod,nsamps);
    dat_Sp = dat_Sp(:,[3,1,2]); % in Z,R,T
    tt01 = mean(tt_s(dat_Sp(:,1)==min(dat_Sp(:,1))));% estimate main P-wave arrival time from first big downswing
    tt_s = tt_s - tt01; tt_s = round_level(tt_s,0.001);
    inwind = (tt_s >= wdo(1)) & (tt_s < wdo(2)); 
    dat_Sp = dat_Sp(inwind,:);
    tt_s = tt_s(inwind);
    try_dat_Sp = flat_hanning_win(tt_s,dat_Sp,wdo(1),wdo(2),3); % 3s taper
    try_dat_Sp = try_dat_Sp/max(max(abs(try_dat_Sp)));
    try_nrm_ZR_sp = norm(try_dat_Sp(:,1)) + norm(try_dat_Sp(:,2));
    % transform to P-SV
    [try_dat_Ps_PSV(:,1),try_dat_Ps_PSV(:,2)] = Rotate_XZ_to_PSV(dat_Ps(:,2),-dat_Ps(:,1),Vp(1),Vs(1),raypp);
    [try_dat_Sp_PSV(:,1),try_dat_Sp_PSV(:,2)] = Rotate_XZ_to_PSV(dat_Sp(:,2),-dat_Sp(:,1),Vp(1),Vs(1),rayps);
    try_dat_Ps_PSV = try_dat_Ps_PSV./max(max(abs(try_dat_Ps_PSV)));
    try_dat_Sp_PSV = try_dat_Sp_PSV./max(max(abs(try_dat_Sp_PSV)));
    
    try_nrm_PSV_ps = norm(try_dat_Ps_PSV(:,1)) + norm(try_dat_Ps_PSV(:,2));
    try_nrm_PSV_sp = norm(try_dat_Sp_PSV(:,1)) + norm(try_dat_Sp_PSV(:,2));

    %% TWO MISFITS
    % ZR
    misfit_ZR(imoh,ivc,1) = xconv_misfit(true_dat_Ps(:,1),true_dat_Ps(:,2),...
                                     try_dat_Ps(:,1) ,try_dat_Ps(:,2))/...
                                     tru_nrm_ZR_ps/try_nrm_ZR_ps;
    misfit_ZR(imoh,ivc,2) = xconv_misfit(true_dat_Sp(:,1),true_dat_Sp(:,2),...
                                     try_dat_Sp(:,1) ,try_dat_Sp(:,2))/...
                                     tru_nrm_ZR_sp/try_nrm_ZR_sp;
    % PSV
    misfit_PSV(imoh,ivc,1) = xconv_misfit(true_dat_Ps_PSV(:,1),true_dat_Ps_PSV(:,2),...
                                     try_dat_Ps_PSV(:,1) ,try_dat_Ps_PSV(:,2))/...
                                     tru_nrm_PSV_ps/try_nrm_PSV_ps;
    misfit_PSV(imoh,ivc,2) = xconv_misfit(true_dat_Sp_PSV(:,1),true_dat_Sp_PSV(:,2),...
                                     try_dat_Sp_PSV(:,1) ,try_dat_Sp_PSV(:,2))/...
                                     tru_nrm_PSV_sp/try_nrm_PSV_sp;

end
end
figure(1), clf
subplot(3,2,1)
contourf(vpvs_try,moh_try,misfit_ZR(:,:,1),30,'edgecolor','none'), colorbar
[~,x,y] = mingrid(misfit_ZR(:,:,1)); hold on, scatter(vpvs_try(x),moh_try(y),70,'w','filled')
ylabel('Zmoh'), title('P-misfit-ZR')

subplot(3,2,3)
contourf(vpvs_try,moh_try,misfit_ZR(:,:,2),30,'edgecolor','none')
[~,x,y] = mingrid(misfit_ZR(:,:,2)); hold on, scatter(vpvs_try(x),moh_try(y),70,'w','filled')
ylabel('Zmoh'), colorbar, title('S-misfit-ZR') 

subplot(3,2,5) 
contourf(vpvs_try,moh_try,sum(misfit_ZR,3),30,'edgecolor','none'), colorbar
[~,x,y] = mingrid(sum(misfit_ZR,3)); hold on, scatter(vpvs_try(x),moh_try(y),70,'w','filled')
xlabel('Vp/Vs'),ylabel('Zmoh'), title('Both-misfit-ZR')

subplot(3,2,2)
contourf(vpvs_try,moh_try,misfit_PSV(:,:,1),30,'edgecolor','none'), colorbar
[~,x,y] = mingrid(misfit_PSV(:,:,1)); hold on, scatter(vpvs_try(x),moh_try(y),70,'w','filled')
title('P-misfit-PSV') 

subplot(3,2,4)
contourf(vpvs_try,moh_try,misfit_PSV(:,:,2),30,'edgecolor','none'), colorbar
[~,x,y] = mingrid(misfit_PSV(:,:,2)); hold on, scatter(vpvs_try(x),moh_try(y),70,'w','filled')
title('S-misfit-PSV') 

subplot(3,2,6)
contourf(vpvs_try,moh_try,sum(misfit_PSV,3),30,'edgecolor','none'), colorbar
[~,x,y] = mingrid(sum(misfit_PSV,3)); hold on, scatter(vpvs_try(x),moh_try(y),70,'w','filled')
xlabel('Vp/Vs'), title('Both-misfit-PSV') 

figure(2);
subplot(2,2,1), plot(tt_p,true_dat_Ps(:,1:2)), ylabel('Ps - ZR')
subplot(2,2,2), plot(tt_p,true_dat_Ps_PSV(:,1:2)), ylabel('Ps - PSV')
subplot(2,2,3), plot(tt_s,true_dat_Sp(:,1:2)), ylabel('Sp - ZR')
subplot(2,2,4), plot(tt_s,true_dat_Sp_PSV(:,1:2)), ylabel('Sp - PSV')

return

figure(3),clf, hold on
plot(moh_try,sum(misfit_ZR,3))
plot(moh_try,sum(misfit_PSV,3))
plot(moh,min(sum(misfit_PSV,3)),'ok')
legend('ZR','PSV')
xlabel('moho depth'), ylabel('misfit')




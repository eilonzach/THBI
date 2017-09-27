% script to test the minimising P-SV correlation technique of estimating
% surface velocity

%% parms
inc = 20;
ph = 'Ps';
synthperiod = 3;
samprate = 20;
nsamps = 2^12;
wdo = [-30 30];
P_err_upwt = 1;

%% model
moh = 42;
Z  = [linspace(0,moh,10) linspace(moh,250,5) 260]';
Vs = [linspace(3.3,3.4,10) linspace(3.9,4.2,5) 4.2]';
Vp = 1.8*Vs;
rho = mantle_vs2rho(Vs,Z);
[zlayt,zlayb,Vs_lay,Vp_lay,rho_lay] = layerise(Z,Vs,0.01,0,Vp,rho);
laymod = struct('nlay',length(Vs_lay),'zlayt',zlayt,'zlayb',zlayb,'Vp',Vp_lay,'Vs',Vs_lay,'rho',rho_lay);
raypp = sind(inc)/Vp(end);
rayps = sind(inc)/Vs(end);


%% Ps synth
[dat_Ps,tt_p] = run_propmat(laymod,'test1','Ps',samprate,inc,synthperiod,nsamps);
dat_Ps = dat_Ps(:,[3,1,2]); % in Z,R,T
tt01 = mean(tt_p(dat_Ps(:,1)==min(dat_Ps(:,1))));% estimate main P-wave arrival time from first big downswing
tt_p = tt_p - tt01; tt_p = round_level(tt_p,0.001);
inwind = (tt_p >= wdo(1)) & (tt_p <= wdo(2)); 
dat_Ps = dat_Ps(inwind,:);
tt_p = tt_p(inwind);
dat_Ps = flat_hanning_win(tt_p,dat_Ps,wdo(1),wdo(2),3); % 3s taper
pmax = max(max(abs(dat_Ps)));
% transform to P-SV
[dat_Ps_PSV(:,1),dat_Ps_PSV(:,2)] = Rotate_XZ_to_PSV(dat_Ps(:,2),-dat_Ps(:,1),Vp(1),Vs(1),raypp);



%% Sp synth
[dat_Sp,tt_s] = run_propmat(laymod,'test1','Sp',samprate,inc,synthperiod,nsamps);

dat_Sp = dat_Sp(:,[3,1,2]); % in Z,R,T
tt01 = mean(tt_s(dat_Sp(:,1)==min(dat_Sp(:,1))));% estimate main P-wave arrival time from first big downswing
tt_s = tt_s - tt01; tt_s = round_level(tt_s,0.001);
inwind = (tt_s >= wdo(1)) & (tt_s <= wdo(2)); 
dat_Sp = dat_Sp(inwind,:);
tt_s = tt_s(inwind);
dat_Sp = flat_hanning_win(tt_s,dat_Sp,wdo(1),wdo(2),3); % 3s taper
smax = max(max(abs(dat_Sp)));
% transform to P-SV
[dat_Sp_PSV(:,1),dat_Sp_PSV(:,2)] = Rotate_XZ_to_PSV(dat_Sp(:,2),-dat_Sp(:,1),Vp(1),Vs(1),rayps);


%% Estimate velocities
[ vp_est_ps,vs_est_ps,E_surf_ps ] = surf_vp_vs_est( dat_Ps,tt_p,raypp,[-5 5],1);
[ vp_est_sp,vs_est_sp,E_surf_sp ] = surf_vp_vs_est( dat_Sp,tt_s,rayps,[-5 5],1);
% combine
vp_range = [2. : 0.025 : 7.5];
vs_range = [1.5 : 0.04 : 4.7];
E_surf_both = E_surf_ps./mingrid(E_surf_ps) + E_surf_sp./mingrid(E_surf_sp)/P_err_upwt;
[~,x,y] = mingrid(E_surf_both);
vs_est_both = vs_range(x);
vp_est_both = vp_range(y);

[dat_Ps_PSV_est(:,1),dat_Ps_PSV_est(:,2)] = Rotate_XZ_to_PSV(dat_Ps(:,2),-dat_Ps(:,1),vp_est_both,vs_est_both,raypp);
[dat_Sp_PSV_est(:,1),dat_Sp_PSV_est(:,2)] = Rotate_XZ_to_PSV(dat_Sp(:,2),-dat_Sp(:,1),vp_est_both,vs_est_both,rayps);



%% plots
figure(11)
% true Vertical and Horizontals
subplot(4,2,1), plot(tt_p,dat_Ps(:,1),'k','linewidth',2), ylabel('Z-ps'), ylim(pmax*[-1 1])
subplot(4,2,3), plot(tt_p,dat_Ps(:,2),'r','linewidth',2), ylabel('R-ps'), ylim(pmax*[-1 1])
subplot(4,2,2), plot(tt_s,dat_Sp(:,1),'k','linewidth',2), ylabel('Z-sp'), ylim(smax*[-1 1])
subplot(4,2,4), plot(tt_s,dat_Sp(:,2),'r','linewidth',2), ylabel('R-sp'), ylim(smax*[-1 1])
% true P and SV
subplot(4,2,5), plot(tt_p,dat_Ps_PSV(:,1),'k','linewidth',2.5), ylabel('P-ps'),  ylim(0.5*pmax*[-1 1])
subplot(4,2,7), plot(tt_p,dat_Ps_PSV(:,2),'r','linewidth',2.5), ylabel('SV-ps'), ylim(0.5*pmax*[-1 1])
subplot(4,2,6), plot(tt_s,dat_Sp_PSV(:,1),'k','linewidth',2.5), ylabel('P-sp'),  ylim(0.5*smax*[-1 1])
subplot(4,2,8), plot(tt_s,dat_Sp_PSV(:,2),'r','linewidth',2.5), ylabel('SV-sp'), ylim(0.5*smax*[-1 1])
% est P and SV
subplot(4,2,5), hold on, plot(tt_p,dat_Ps_PSV_est(:,1),'--c')
subplot(4,2,7), hold on, plot(tt_p,dat_Ps_PSV_est(:,2),'--g')
subplot(4,2,6), hold on, plot(tt_s,dat_Sp_PSV_est(:,1),'--c')
subplot(4,2,8), hold on, plot(tt_s,dat_Sp_PSV_est(:,2),'--g')

figure(12)
subplot(1,3,1), hold on, title('Ps error map')
    contourf(vs_range,vp_range,E_surf_ps,40,'edgecolor','none'), colorbar
    scatter(vs_est_ps,vp_est_ps,50,'r','filled')
subplot(1,3,2), hold on, title('Sp error map')
    contourf(vs_range,vp_range,E_surf_sp,40,'edgecolor','none'), colorbar
    scatter(vs_est_sp,vp_est_sp,50,'g','filled')
subplot(1,3,3), hold on, title('Combined error map')
    contourf(vs_range,vp_range,E_surf_both,40,'edgecolor','none'), colorbar
    scatter(vs_est_ps,vp_est_ps,70,'r','linewidth',2)
    scatter(vs_est_sp,vp_est_sp,70,'g','linewidth',2)
    scatter(Vs(1),Vp(1),70,'w','linewidth',2)
    scatter(vs_est_both,vp_est_both,50,'w','filled')






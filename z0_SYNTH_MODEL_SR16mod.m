function z0_SYNTH_MODEL_SR16mod(par)

global TRUEmodel TLM


%% GET TARGET MODEL
plat = 42;
plon = 254;
% read 1D profile from Shen and Ritzwoller model
ncfile = '~/Work/data/models_seismic/SR16_3d_Vs/US.2016.nc';
[ TRUEmodel ] = read_seismodel_nc( ncfile );

TRUEmodel.vs = squeeze(TRUEmodel.Vsv(mindex(TRUEmodel.lat,plat),mindex(TRUEmodel.lon,mod(plon,360)),:));
TRUEmodel.vp = squeeze(TRUEmodel.Vpv(mindex(TRUEmodel.lat,plat),mindex(TRUEmodel.lon,mod(plon,360)),:));
TRUEmodel.rho = squeeze(TRUEmodel.rho(mindex(TRUEmodel.lat,plat),mindex(TRUEmodel.lon,mod(plon,360)),:));

% % RHO is effed up = we don't agree on rho scaling, but it matters a lot for
% % the elastic moduli that propmat actually uses. 
% rho = linterp(trymod.z,trymod.rho,Z);
% % vs = linterp(trymod.z,trymod.VS,Z)+0.3;
% vp = linterp(trymod.z,trymod.VP,Z);


%% MAKE LAYERISED TARGET MODEL
% layerise
[zlayt,zlayb,Vslay,Vplay,rholay] = ...
    layerise(TRUEmodel.Z,TRUEmodel.vs,par.forc.mindV,0,TRUEmodel.vp,TRUEmodel.rho); 
nlay = length(Vslay);
TRUElaymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay);
TLM = TRUElaymodel;


%% PLOT FINAL MODEL
figure(95); clf; set(gcf,'pos',[120 151 920 947])

subplot(131), hold on;
plot(TRUEmodel.Vp,TRUEmodel.Z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);

subplot(132), hold on;
plot(TRUEmodel.Vs,TRUEmodel.Z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);

subplot(133), hold on;
plot(TRUEmodel.rho,TRUEmodel.Z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);



end


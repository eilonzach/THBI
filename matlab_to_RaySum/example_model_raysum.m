% clear
% close

mindV = 0.08;
% pick a location
plat = 43;
plon = 247;



profile clear
profile on
addpath('../functions');
% read 1D profile from Shen and Ritzwoller model
ncfile = '~/Work/data/models_seismic/SR16_3d_Vs/US.2016.nc';
[ model ] = read_seismodel_nc( ncfile );

vs = squeeze(model.Vsv(mindex(model.lat,plat),mindex(model.lon,mod(plon,360)),:));
vp = squeeze(model.Vpv(mindex(model.lat,plat),mindex(model.lon,mod(plon,360)),:));
rho = squeeze(model.rho(mindex(model.lat,plat),mindex(model.lon,mod(plon,360)),:));
Z = model.Z;

% layerise profile
% [zlayt,zlayb,Vslay,Vplay,rho] = layerise(Z,vs,mindV,1,vp,rho); nlay = length(Vslay);
[zlayt,zlayb,Vslay] = layerise(Z,vs,mindV,0); nlay = length(Vslay);
if zlayb(end)==max(Z); zlayb(end)=zlayt(end); end % half space
Vplay = 1.81*Vslay;
rho = sqrt(Vplay);
fprintf('Nlay = %.0f\n',nlay);

%% plot
figure(1); clf, hold on
plot(vs,Z,'-ko')
plot(vp,Z,'-ko')
zlayp = reshape([zlayt';zlayb'],2*nlay,1);
vslayp = reshape([Vslay';Vslay'],2*nlay,1);
vplayp = reshape([Vplay';Vplay'],2*nlay,1);
plot(vslayp,zlayp,'-ro')
plot(vplayp,zlayp,'-ro')
set(gca,'ydir','reverse','ylim',[0, max(Z)],'xlim',[0.9*min(vs) 1.1*max(vp)])


model1D = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rho);

% write to Raysum format
writeRAYSUMmodfile( model1D,'EGrsinputfile.mod')
writeRAYSUMgeomfile( [0:30:90]',[2.5:0.5:4]'*1e-5,'EGrsinputfile.geom' )
% do raysum on it
system('~/Work/codes/Raysum_fwd_v1.2/bin/seis-spread EGrsinputfile.mod EGrsinputfile.geom out.ph out.arr out.tr')
% read raysum output
[traces,tt] = readraysumtr('out.tr');
% plot
figure(2); clf, hold on
for ip = 1:3
subplot(3,1,ip)
plot(tt,squeeze(traces(:,ip,:)),'Linewidth',1.5)
end
% profile viewer
function [model,laymodel] = z0_SYNTH_MODEL_LAB_TEST(par,zsed,zmoh,zlab,wlab,flab,ifplot)
% [model,laymodel] = z0_SYNTH_MODEL_LAB_TEST(par,zsed,zmoh,zlab,wlab,flab,ifplot)


if nargin < 7 || isempty(ifplot)
    ifplot=false;
end

global TRUEmodel TLM

dvdz_cr = 0.008; % Vs increase with depth in the crust, (km/s)/km
dvdz_lab = 0.003; % Vs increase with depth in the mantle, (km/s)/km
dvdz_asth = 0.001; % Vs increase with depth in the mantle, (km/s)/km

% parms
% zsed = 0;
% zmoh = 50;
vpvs = 1.8;
mxi = 1.04;
cxi = 0.95;
zlabm = zlab;
% wlab = 10;
dz = par.mod.dz;
zmax = par.mod.maxz;
% flab = 0.1;

%mantle velocities
vlt = 4.3; % top of lith
vlab = 4.3 + dvdz_lab*(zlabm-zmoh); % mean v at lab
vlb = vlab*(1 + 0.5*flab); % v at top of lab
vat = vlab*(1 - 0.5*flab); % v at bot of lab
vab = vat + dvdz_asth*(zmax-zmoh); % v at bot of lab


% flat seds
zs =  [ 0 zsed]'; if zsed == 0, zs = []; end
vss = [2.9 3.1]'; if zsed == 0, vss= []; end
% ramp crust
zc =  unique([zsed:dz:zmoh,zmoh])';
vsc = interp1([zsed,zmoh],3.5+[0 dvdz_cr*(zmoh-zsed)]',zc);
% shallow ramp up in mantle lith
zl =  unique([zmoh:dz:zlabm-wlab/2,zlabm-wlab/2])';
vsl = interp1([zmoh,zlabm-wlab/2],[vlt,vlb],zl);
% steep ramp down in LAB
if wlab == 0
    zlab = [];
    vslab = [];
else
    zlab =  unique([zlabm-wlab/2:dz:zlabm+wlab/2,zlabm+wlab/2])';
    vslab = interp1(zlabm+wlab*[-0.5 +0.5],[vlb,vat],zlab);
end
% ramp up in asth
za =  unique([zlabm+wlab/2:dz:zmax,zmax])';
vsa = interp1([zlabm+wlab/2,zmax],[vat,vab],za);

% %random knick point
% vsl(zl<90 & zl>80) = 3.9;


Z  = [zs; zc; zl; zlab; za];
Vs = [vss;vsc;vsl;vslab;vsa];
Vp = [sed_vs2vp(vss); vpvs*vsc;mantle_vs2vp([vsl;vslab;vsa],[zl; zlab; za])];
rho= [sed_vs2rho([vss;vsc]);mantle_vs2rho([vsl;vslab;vsa],[zl; zlab; za])];

%% Put into TRUEMODEL struct
TRUEmodel = struct([]);
TRUEmodel(1).z = Z;
TRUEmodel.VS = Vs;
TRUEmodel.VP = Vp;
TRUEmodel.rho = rho;
TRUEmodel.vpvs = vpvs;
TRUEmodel.cxi = cxi;
TRUEmodel.mxi = mxi;
TRUEmodel.Nz = length(Z);
TRUEmodel.zsed = zsed;
TRUEmodel.zmoh = zmoh;
TRUEmodel.selev = 0;
TRUEmodel.Sanis = zeros(length(Z),1);
TRUEmodel.Panis = zeros(length(Z),1);
TRUEmodel.Z = TRUEmodel.z;
TRUEmodel.vs = TRUEmodel.VS;
TRUEmodel.vp = TRUEmodel.VP;

%% ===================  LAYERISE PROFILE  ===================
[zlayt,zlayb,Vslay] = ...
    layerise(Z,Vs,par.forc.mindV/3,0); 
nlay = length(Vslay);

% S to P and rho structure
xs = 1:find(zlayb==TRUEmodel.zsed); if TRUEmodel.zsed ==0, xs = []; end
xc = find(zlayt==TRUEmodel.zsed):find(zlayb==TRUEmodel.zmoh);
xm = find(zlayt==TRUEmodel.zmoh):nlay;
Vplay = [sed_vs2vp(Vslay(xs));...
         TRUEmodel.vpvs*Vslay(xc);...
         mantle_vs2vp(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
rholay = [sed_vs2rho(Vslay([xs,xc]));...
          mantle_vs2rho(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
xilay = [zeros(length(xs),1);...
         TRUEmodel.cxi*ones(length(xc),1);...
         TRUEmodel.mxi*ones(length(xm),1)]; % S radial anisotropy
philay = ones(nlay,1); % P radial anisotropy
etalay = ones(nlay,1); % eta anisotropy

TLM = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay,'xi',xilay,'phi',philay,'eta',etalay);
if any(isnan(TLM.rho))
    error('NaN densities')
end




%% PLOT FINAL MODEL
if ifplot
figure(95); clf; set(gcf,'pos',[120 151 920 947])
subplot(131), hold on;
plot(TRUEmodel.vp,TRUEmodel.z,'-k','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
subplot(132), hold on;
plot(TRUEmodel.vs,TRUEmodel.z,'-k','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
subplot(133), hold on;
plot(TRUEmodel.rho,TRUEmodel.z,'-k','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
end

model = TRUEmodel;
laymodel = TLM;

end


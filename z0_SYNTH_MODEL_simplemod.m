function [model,laymodel] = z0_SYNTH_MODEL_simplemod(par,ifplot)


if nargin < 2 || isempty(ifplot)
    ifplot=false;
end

global TRUEmodel TLM


% parms
hsed = 5;
H = 50;
K = 1.8;
LAB = 125;
LABgradz = 20;


% % simple layered model
% zlayt = [0  15 20  55 110 160 220 260]';
% zlayb = [15 20 55 110 160 220 260 300]';
% Vs_lay = [3.2 3.8 4.5 4.7 4.4 4.7 4.25 4.95]';

% zlayt = [0  H  LAB ]';
% zlayb = [H  LAB 200 ]';
% Vs_lay = [3.7 4.5  4.2 ]';
% Vp_lay = K*Vs_lay;
% rho_lay = mantle_vs2rho(Vs_lay,0.5*(zlayb+zlayt));


% Z  = [  0   H   H  LAB linspace(LAB,200,20)]';
% Vs = [3.7 3.7 4.5  4.5 linspace(4.2,4.5,20)]';

% flat seds
zs =  [ 0 hsed]';
vss = [3.2 3.2]';
% flat crust
zc =  [hsed  H]';
vsc = [3.7 3.7]';
% shallow ramp up in mantle lith
zl =  unique([H:par.mod.dz:LAB-LABgradz/2,LAB-LABgradz/2])';
vsl = interp1([H,LAB-LABgradz/2],[4.3,4.48],zl);
% steep ramp down in LAB
zlab =  unique([LAB-LABgradz/2:par.mod.dz:LAB+LABgradz/2,LAB+LABgradz/2])';
vslab = interp1(LAB+LABgradz*[-0.5 +0.5],[4.48,4.20],zlab);
% ramp up in asth
za =  unique([LAB+LABgradz/2:par.mod.dz:200,200])';
vsa = interp1([LAB+LABgradz/2,200],[4.20,4.5],za);

% %random knick point
% vsl(zl<90 & zl>80) = 3.9;


Z  = [zs; zc; zl; zlab; za];
Vs = [vss;vsc;vsl;vslab;vsa];
Vp = [sed_vs2vp(vss);K*vsc;K*vsl;K*vslab;K*vsa];
rho= [sed_vs2rho([vss;vsc]);mantle_vs2rho([vsl;vslab;vsa],[zl; zlab; za])];

%% Put into TRUEMODEL struct
TRUEmodel = struct([]);
TRUEmodel(1).z = Z;
TRUEmodel.VS = Vs;
TRUEmodel.VP = Vp;
TRUEmodel.rho = rho;
TRUEmodel.Nz = length(Z);
TRUEmodel.zsed = hsed;
TRUEmodel.zmoh = H;
TRUEmodel.selev = 0;
TRUEmodel.Sanis = zeros(length(Z),1);
TRUEmodel.Panis = zeros(length(Z),1);
TRUEmodel.Z = TRUEmodel.z;
TRUEmodel.vs = TRUEmodel.VS;
TRUEmodel.vp = TRUEmodel.VP;

%% MAKE LAYERISED TARGET MODEL
% layerise
[zlayt,zlayb,Vs_lay,Vp_lay,rho_lay] = layerise(Z,Vs,0.05,1,Vp,rho);
TLM = struct('nlay',length(Vs_lay),'zlayt',zlayt,'zlayb',zlayb,'Vp',Vp_lay,'Vs',Vs_lay,'rho',rho_lay);
           




%% PLOT FINAL MODEL
if ifplot
figure(95); clf; set(gcf,'pos',[120 151 920 947])
subplot(131), hold on;
plot(TRUEmodel.vp,TRUEmodel.z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
subplot(132), hold on;
plot(TRUEmodel.vs,TRUEmodel.z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
subplot(133), hold on;
plot(TRUEmodel.rho,TRUEmodel.z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
end

model = TRUEmodel;
laymodel = TLM;

end


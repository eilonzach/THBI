function [model,laymodel] = z0_SYNTH_MODEL_custommod(par,ifplot)
%  [model,laymodel] = z0_SYNTH_MODEL_custommod(par,ifplot)

if nargin < 2 || isempty(ifplot)
    ifplot=false;
end

global TRUEmodel TLM


%% resolve important values from prior distributions
selev =0;
h_sed = 5;
h_crust = 50.5;

k_crust = 4;
k_mantle = 6;

vs_sed = [3.2 3.2];
kvs_crust = [3.5 3.7 3.77 3.81];
kvs_mantle = [4.1 4.4 4.4 4.05 4.2 4.4];

% sharper LAB
k_mantle = 9;
kvs_mantle = [4.2 4.4 4.4 4.45 4.5 4.0 4.1 4.2 4.45];

vpvs_crust = 1.8;

xi_crust = 1.04;
xi_mantle = 1.05;



%% DERIVATIVE PARMS
% DEPTHS
cminz = h_sed;
cmaxz = h_sed+h_crust;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
mminz = h_sed+h_crust;
mmaxz = par.mod.maxz + selev;
zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';
% CRUST splines
cknots = linspace(cminz,cmaxz,k_crust-1)';
[cspbasis] = make_splines(cknots,par,zc,zc);
% MANTLE splines
mknots = linspace(mminz,mmaxz,k_mantle-1)';
[mspbasis] = make_splines(mknots,par,zm,zm);
% OVERALL
M = 1 + 2 + 1 + k_crust + k_mantle + 1 + 2;


%% MAKE ALL PARAMETER STRUCTURES
sed = struct('h',h_sed,'VS',vs_sed);
crust = struct('h',h_crust,'Nsp',k_crust+1,'Nkn',k_crust,'VS_sp',kvs_crust,'vpvs',vpvs_crust,'xi',xi_crust,'splines',cspbasis,'knots',cknots,'z_sp',zc);
mantle = struct('Nsp',k_mantle+1,'Nkn',k_mantle,'VS_sp',kvs_mantle,'xi',xi_mantle,'splines',mspbasis,'knots',mknots,'z_sp',zm);
data = struct('sig_Ps_RF',par.mod.data.prior_sigma.BW.Ps,...
              'sig_Sp_RF',par.mod.data.prior_sigma.BW.Ps,...
              'sig_SW',par.mod.data.prior_sigma.SW);


%% MODEL WITH ALL PARMS
model = struct('sedmparm',sed,'crustmparm',crust,'mantmparm',mantle,...
               'datahparm',data,'M',M,'selev',selev);
           
%% TURN PARMS INTO REAL TARGET MODEL
TRUEmodel = make_mod_from_parms(model,par);
model = TRUEmodel;
% same format...

TRUEmodel.Z = TRUEmodel.z;
TRUEmodel.vs = TRUEmodel.VS;
TRUEmodel.vp = TRUEmodel.VP;

% % RHO is effed up = we don't agree on rho scaling, but it matters a lot for
% % the elastic moduli that propmat actually uses. 
% rho = linterp(trymod.z,trymod.rho,Z);
% % vs = linterp(trymod.z,trymod.VS,Z)+0.3;
% vp = linterp(trymod.z,trymod.VP,Z);


%% MAKE LAYERISED TARGET MODEL
% layerise
[zlayt,zlayb,Vslay,Vplay,rholay] = ...
    layerise(TRUEmodel.z,TRUEmodel.vs,par.forc.mindV,0,TRUEmodel.vp,TRUEmodel.rho); 
nlay = length(Vslay);
laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay);
TLM = laymodel;


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



end


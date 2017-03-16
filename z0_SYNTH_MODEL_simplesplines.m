    function [model,laymodel] = z0_SYNTH_MODEL_simplesplines(par,ifplot)


if nargin < 2 || isempty(ifplot)
    ifplot=false;
end

global TRUEmodel TLM


%% CHOOSE CUSTOM KEY PARAMETERS
selev =0;
h_sed = 5;
h_crust = 50;

k_crust = 4;

vs_sed = [3.2 3.2];
kvs_crust = [3.5 3.7 3.77 3.81];
kvs_mantle = [4.1 4.4 4.4 4.05 4.2 4.4];

% sharper LAB
k_mantle = 6;
kvs_mantle = [4.2 4.4  4.5 4.0 4.1 4.45];

vpvs_crust = 1.8;


%% DERIVATIVE PARMS
% DEPTHS
cminz = h_sed;
cmaxz = h_sed+h_crust;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
mminz = h_sed+h_crust;
mmaxz = par.mod.maxz + selev;
zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';
% CRUST splines
dzsp = (cmaxz-cminz)/(k_crust-2);
knots = [repmat(cminz,1,3),cminz:dzsp:cmaxz,repmat(cmaxz,1,3)]';
sp = fastBSpline.lsqspline(knots,2,linspace(cminz,cmaxz,k_crust)',kvs_crust); % dummy velocities as placeholder
cspbasis = sp.getBasis(zc); cspbasis = cspbasis(:,2:end-1);
% MANTLE splines
dzsp = (mmaxz-mminz)/(k_mantle-2);
knots = [repmat(mminz,1,3),mminz:dzsp:mmaxz,repmat(mmaxz,1,3)]';
sp = fastBSpline.lsqspline(knots,2,linspace(mminz,mmaxz,k_mantle)',kvs_mantle); % dummy velocities as placeholder
mspbasis = sp.getBasis(zm); mspbasis = mspbasis(:,2:end-1);
% OVERALL
M = 1 + 2 + 1 + k_crust + k_mantle + 1;


%% MAKE ALL PARAMETER STRUCTURES
sed = struct('h',h_sed,'VS',vs_sed);
crust = struct('h',h_crust,'Nkn',k_crust,'VS_sp',kvs_crust,'splines',cspbasis,'vpvs',vpvs_crust);
mantle = struct('Nkn',k_mantle,'VS_sp',kvs_mantle,'splines',mspbasis);
data = struct('sigmaPsRF',par.mod.data.prior_sigmaPsRF,...
              'sigmaSpRF',par.mod.data.prior_sigmaSpRF,...
              'sigmaSW',par.mod.data.prior_sigmaSW);


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


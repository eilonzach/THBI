function [ final_model ] = c4_FINAL_MODEL( posterior,allmodels,par,ifsave,ofile )
% [ final_model ] = c4_FINAL_MODEL( posterior,allmodels,par,ifsave )
%   Detailed explanation goes here

if nargin < 4 || isempty(ifsave)
    ifsave = 0;
end

if nargin < 5 || isempty(ofile)
    ofile = './final_model';
end

o4 = ones(1,4);

% don't use pre burn-in models!
allmodels([allmodels.bestmods]'==false) =  [];

fprintf('  > Gathering models (layer-wise)\n    & fitting gaussians to each parm/depth\n')

%% ===================== GATHER MODELS, FIT WITH GAUSSIANS ================
%% DISCONTINUITIES
Zd = struct('mu',[],'std',[]);
% crust/sed
X = midpts([0:0.1:par.mod.sed.hmax]);
Nds = hist(posterior.zsed,X);
if sum(Nds~=0)==1
    Zd(1).mu = mean(posterior.zsed); Zd(1).std = 0.001;
else
    [Zd(1).std,Zd(1).mu] = gaussfit( X, Nds );
end

% Moho
X = midpts([0:0.5:par.mod.sed.hmax+par.mod.crust.hmax]);
Nds = hist(posterior.zmoh,X);
if sum(Nds~=0)==1
    Zd(2).mu = mean(posterior.zmoh); Zd(2).std = 0.001;
else
    [Zd(2).std,Zd(2).mu] = gaussfit( X, Nds );
end   
%% SEDIMENTS
Zs = [0;Zd(1).mu];
VSs = struct('mu',zeros(2,1),'std',zeros(2,1));   
X = midpts([par.mod.sed.vsmin:0.05:par.mod.sed.vsmax]);
Lvss1 = hist(posterior.VSsedtop,X);
Lvss2 = hist(posterior.VSsedbot,X);
if sum(Lvss1~=0)==1
    VSs.mu(1) = mean(posterior.VSsedtop); VSs.std(1) = 0.001;
else
    [VSs.std(1),VSs.mu(1)] = gaussfit(X,Lvss1);
end
if sum(Lvss2~=0)==1
    VSs.mu(2) = mean(posterior.VSsedbot); VSs.std(2) = 0.001;
else
    [VSs.std(2),VSs.mu(2)] = gaussfit(X,Lvss2);
end
% Vp
VPs.mu = sed_vs2vp(VSs.mu);
VPs.std = 0.5*[sed_vs2vp(VSs.mu+VSs.std)-sed_vs2vp(VSs.mu-VSs.std)];
%rho
rhos.mu  = sed_vs2rho(VSs.mu);
rhos.std = 0.5*[sed_vs2rho(VSs.mu+VSs.std)-sed_vs2rho(VSs.mu-VSs.std)];

%% VPVS
X = midpts([par.mod.crust.vpvsmin:0.01:par.mod.crust.vpvsmax]);
Lvpvs = hist(posterior.vpvs,X);
vpvs = struct('mu',zeros(2,1),'std',zeros(2,1));   
if sum(Lvpvs~=0)==1
    vpvs.mu = mean(posterior.vpvs); vpvs.std = NaN;
else
    [vpvs.std,vpvs.mu] = gaussfit(X,Lvpvs);
end

%% CRUST
% top and bottom
X = midpts([par.mod.crust.vsmin:0.05:par.mod.crust.vsmax]);
XX = [par.mod.crust.vsmin*o4, X ,par.mod.crust.vsmax*o4];

VSctb = struct('mu',zeros(2,1),'std',zeros(2,1));   
Lvsct = hist(posterior.VScrusttop,X); LLvsct = [Lvsct(1)*o4,Lvsct,Lvsct(end)*o4];
Lvscb = hist(posterior.VScrustbot,X); LLvscb = [Lvscb(1)*o4,Lvscb,Lvscb(end)*o4];
if sum(Lvsct~=0)==1
    VSctb.mu(1) = mean(posterior.VScrusttop); VSctb.std(1) = 0.001;
else
    [VSctb.std(1),VSctb.mu(1)] = gaussfit( XX, LLvsct );
end
if sum(Lvscb~=0)==1
    VSctb.mu(2) = mean(posterior.VScrustbot); VSctb.std(2) = 0.001;
else
    [VSctb.std(2),VSctb.mu(2)] = gaussfit(XX,LLvscb);
end
% middle points (more than 1std from each edge
% need nansum in case the edge discontinuities have zero std
Zcg = unique([nansum([Zd(1).mu,Zd(1).std]): 1 : nansum([Zd(2).mu,-Zd(2).std])]');
Nzcg = length(Zcg);
% gather velocities interpolated on these points
VScg = zeros(Nzcg,posterior.Nstored);
VPcg = zeros(Nzcg,posterior.Nstored);
for ii = 1:posterior.Nstored
    VScg(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VS,Zcg);
    VPcg(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VP,Zcg);
end
% fit gaussians at each depth of the gather
VScz = struct('mu',zeros(Nzcg,1),'std',zeros(Nzcg,1));
VPcz = struct('mu',zeros(Nzcg,1),'std',zeros(Nzcg,1));
for iz = 1:Nzcg
    Lvsc = hist(VScg(iz,:),X); LLvsc = [Lvsc(1)*o4,Lvsc,Lvsc(end)*o4];
    if sum(Lvsc~=0)==1
        VScz.mu(iz) = mean(VScg(iz,:)); VScz.std(iz) = 0.001;
    else
        [VScz.std(iz), VScz.mu(iz)] = gaussfit( XX, LLvsc );
    end
    
    Lvpc = hist(VPcg(iz,:),vpvs.mu*X); LLvpc = [Lvpc(1)*o4,Lvpc,Lvpc(end)*o4];
    if sum(Lvpc~=0)==1
        VPcz.mu(iz) = mean(VPcg(iz,:)); VPcz.std(iz)=0.001;
    else
        [VPcz.std(iz), VPcz.mu(iz)] = gaussfit( vpvs.mu*XX, LLvpc );
    end
end
% STITCH middle to edges
    % at the layer edges, the mean values at that depth might be
    % susceptible to sampling of the layers on either side - i.e. sampling
    % across the discontinuity. To avoid this, pin the edges to the
    % top/bottom values, and linearly grade towards these values over a
    % length scale defined by the standard deviation of the discontinuity
    % depth. Use erf function for this.
Zc = unique([Zd(1).mu:1:Zd(2).mu,Zd(2).mu]');
Nzc = length(Zc);
VSc = struct('mu',zeros(Nzc,1),'std',zeros(Nzc,1));
VPc = struct('mu',zeros(Nzc,1),'std',zeros(Nzc,1));
% Vs
VSc.mu  = linterp([Zd(1).mu;Zcg;Zd(2).mu],...
                  [VSctb.mu(1);VScz.mu;VSctb.mu(2)],...
                   Zc);
VSc.std = linterp([Zd(1).mu;Zcg;Zd(2).mu],...
                  [VSctb.std(1);VScz.std;VSctb.std(2)],...
                   Zc);
% Vp
VPc.mu  = linterp([Zd(1).mu;Zcg;Zd(2).mu],...
                  [vpvs.mu*VSctb.mu(1);VPcz.mu;vpvs.mu*VSctb.mu(2)],...
                   Zc);
VPc.std = linterp([Zd(1).mu;Zcg;Zd(2).mu],...
                  [vpvs.mu*VSctb.std(1);VPcz.std;vpvs.mu*VSctb.std(2)],...
                   Zc);
% rho
rhoc.mu  = sed_vs2rho(VSc.mu);
rhoc.std = 0.5*[sed_vs2rho(VSc.mu+VSc.std)-sed_vs2rho(VSc.mu-VSc.std)];

%% MANTLE 
% top
X = midpts([par.mod.crust.vsmin:0.05:par.mod.mantle.vsmax]);
XX = [par.mod.crust.vsmin*o4, X ,par.mod.mantle.vsmax*o4];

VSmt = struct('mu',0,'std',0);   
Lvsmt = hist(posterior.VSmanttop,X); LLvsmt = [Lvsmt(1)*o4, Lvsmt, Lvsmt(end)*o4];
if sum(Lvsmt~=0)==1
    VSmt.mu = mean(posterior.VSmanttop); VSmt.std = 0.001;
else
    [VSmt.std,VSmt.mu] = gaussfit(XX,LLvsmt);
end

% middle points (more than 1std from top edge
Zmg = unique([nansum([Zd(2).mu,Zd(2).std]): 2 : par.mod.maxz,par.mod.maxz]');
Nzmg = length(Zmg);

% gather velocities interpolated on these points
VSmg = zeros(Nzmg,posterior.Nstored);
VPmg = zeros(Nzmg,posterior.Nstored);
for ii = 1:posterior.Nstored
    VSmg(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VS,Zmg);
    VPmg(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VP,Zmg);
end
% fit gaussians at each depth of the gather
VSmz = struct('mu',zeros(Nzmg,1),'std',zeros(Nzmg,1));
% VPmz = struct('mu',zeros(Nzmg,1),'std',zeros(Nzmg,1));
for iz = 1:Nzmg
    Lvsm = hist(VSmg(iz,:),X); LLvsm = [Lvsm(1)*o4,Lvsm,Lvsm(end)*o4];
    if sum(Lvsm~=0)==1
        VSmz.mu(iz) = mean(VSmg(iz,:)); VSmz.std(iz) = 0.001;
    else
        [VSmz.std(iz), VSmz.mu(iz)] = gaussfit( XX, LLvsm );
    end
    
%     Lvpm = hist(VPmg(iz,:),mantle_vs2vp(X,Zmg(iz))); LLvpm = [Lvpm(1)*o4,Lvpm,Lvpm(end)*o4];
%     if sum(Lvpm~=0)==1
%         VPmz.mu(iz) = mean(VPmg(iz,:)); VPmz.std(iz)=NaN;
%     else
%         [VPmz.std(iz), VPmz.mu(iz)] = gaussfit( mantle_vs2vp(XX,Zmg(iz)), LLvpm );
%     end
end
% since now Vp is 1:1 function of Vs in the mantle:
VPmz = struct('mu',mantle_vs2vp(VSmz.mu,Zmg),'std',mantle_vs2vp(VSmz.std,Zmg));

% STITCH middle to edges
Zm = unique([Zd(2).mu : 2 : par.mod.maxz,par.mod.maxz]');
Nzm = length(Zm);
VSm = struct('mu',zeros(Nzm,1),'std',zeros(Nzm,1));
VPm = struct('mu',zeros(Nzm,1),'std',zeros(Nzm,1));
% Vs
VSm.mu  = interp1([Zd(2).mu;Zmg],[VSmt.mu; VSmz.mu], Zm);
VSm.std = interp1([Zd(2).mu;Zmg],[VSmt.std;VSmz.std],Zm);
% Vp
VPm.mu  = interp1([Zd(2).mu;Zmg],[mantle_vs2vp(VSmt.mu,Zd(2).mu) ;VPmz.mu], Zm);
VPm.std = interp1([Zd(2).mu;Zmg],[mantle_vs2vp(VSmt.std,Zd(2).mu);VPmz.std],Zm);

% rho
rhom.mu  = mantle_vs2rho(VSm.mu,Zm);
rhom.std = 0.5*[mantle_vs2rho(VSm.mu+VSm.std,Zm)-mantle_vs2rho(VSm.mu-VSm.std,Zm)];

%%  HYPERPARAMETERS 
% sigmaPsRF
postvar = lognfit(posterior.datahparm.sigmaPsRF,0.05);
sigmaPsRF(1).mu_log10 = postvar(1)/log(10);
sigmaPsRF(1).std_log10 = postvar(2)/log(10);
% sigmaSpRF
postvar = lognfit(posterior.datahparm.sigmaSpRF,0.05);
sigmaSpRF(1).mu_log10 = postvar(1)/log(10);
sigmaSpRF(1).std_log10 = postvar(2)/log(10);
% sigmaSW
postvar = lognfit(posterior.datahparm.sigmaSW,0.05);
sigmaSW(1).mu_log10 = postvar(1)/log(10);
sigmaSW(1).std_log10 = postvar(2)/log(10);
% sigmaPsRF_lo
postvar = lognfit(posterior.datahparm.sigmaPsRF_lo,0.05);
sigmaPsRF_lo(1).mu_log10 = postvar(1)/log(10);
sigmaPsRF_lo(1).std_log10 = postvar(2)/log(10);
% sigmaSpRF_lo
postvar = lognfit(posterior.datahparm.sigmaSpRF_lo,0.05);
sigmaSpRF_lo(1).mu_log10 = postvar(1)/log(10);
sigmaSpRF_lo(1).std_log10 = postvar(2)/log(10);

%% ===================== ADD ALL LAYERS TOGETHER ================

final_model.Z = [Zs;Zc;Zm];
final_model.VSbest = [VSs.mu;VSc.mu;VSm.mu];
final_model.VSstd = [VSs.std;VSc.std;VSm.std];
final_model.VPbest = [VPs.mu;VPc.mu;VPm.mu];
final_model.VPstd = [VPs.std;VPc.std;VPm.std];
final_model.rhobest = [rhos.mu;rhoc.mu;rhom.mu];
final_model.rhostd = [rhos.std;rhoc.std;rhom.std];
final_model.hyperparms = struct('sigmaPsRF',sigmaPsRF,'sigmaSpRF',sigmaSpRF,'sigmaSW',sigmaSW,'sigmaPsRF_lo',sigmaPsRF_lo,'sigmaSpRF_lo',sigmaSpRF_lo);

%% =====================  Suite of models  =====================
Zsuite = sort(unique([[0:2:par.mod.maxz]';[-5:0.2:5]'+Zd(1).mu;[-5:0.2:5]'+Zd(2).mu]));
Zsuite(Zsuite<0) = [];
Nz = length(Zsuite);

fprintf('  > Resolving suite of models onto common basis \n')
VSsuite = zeros(Nz,posterior.Nstored);
VPsuite = zeros(Nz,posterior.Nstored);
rhosuite = zeros(Nz,posterior.Nstored);
for ii = 1:posterior.Nstored
VSsuite(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VS,Zsuite);
VPsuite(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VP,Zsuite);
rhosuite(:,ii) = linterp(allmodels(ii).z,allmodels(ii).rho,Zsuite);
end

suite_of_models.Z = Zsuite;
suite_of_models.VS = VSsuite;
suite_of_models.VP = VPsuite;
suite_of_models.rho = rhosuite;



if ifsave
    save(ofile,'final_model')
end




end









function prior = a2_BUILD_PRIOR(par,Niter,zatdep)
% prior = a2_BUILD_PRIOR(par,Niter,zatdep)
% profile clear
% profile on

if nargin <2 || isempty(Niter) 
    Niter = 1e4;
end

o = [];

prior = struct('Niter',Niter','Npass',0,...
               'zsed',o,'zmoh',o,...
               'kcrust',o,'kmantle',o,...
               'VSsedtop',o,'VSsedbot',o,'VScrusttop',o,'VScrustbot',o,...
               'VSmantle',repmat(o,1,6),'zatdep',zeros(6,1),...
               'Zkn_crust',nan(1,par.mod.crust.kmax),'Zkn_mantle',nan(1,par.mod.mantle.kmax),...
               'fdVSsed',o,'fdVSmoh',o,'vpvs',o,'cxi',o,'mxi',o);

if nargin<3 || isempty(zatdep)
    prior.zatdep = linspace(par.mod.sed.hmax+par.mod.crust.hmin+0.1,par.mod.maxz,50)';
else
    prior.zatdep = zatdep;
end

% prep for par
parprior = prior;           
parprior(Niter) = prior(1); % set last val to give prior right dimension

passed = zeros(Niter,1);
           
for kk = 1:Niter
    model = b1_INITIATE_MODEL(par,0,0);
    ifpass = a1_TEST_CONDITIONS( model, par );
    
    if ~ifpass, continue; end
    
    passed(kk) = 1;
    
    parprior(kk).zsed = model.zsed; %#ok<PFOUS>
    parprior(kk).zmoh = model.zmoh;
    parprior(kk).kcrust = model.crustmparm.Nkn;
    parprior(kk).kmantle = model.mantmparm.Nkn;
    parprior(kk).VSsedtop = model.sedmparm.VS(1);
    parprior(kk).VSsedbot = model.sedmparm.VS(2);
    parprior(kk).VScrusttop = model.crustmparm.VS_sp(1);
    parprior(kk).VScrustbot = model.crustmparm.VS_sp(end);
    parprior(kk).VSmantle(1,:) = linterp(model.z,model.VS,prior.zatdep);
    parprior(kk).fdVSsed = model.fdVSsed;
    parprior(kk).fdVSmoh = model.fdVSmoh;
    parprior(kk).Zkn_crust(1:model.crustmparm.Nkn-2) = model.crustmparm.knots(2:end-1);
    parprior(kk).Zkn_mantle(1:model.mantmparm.Nkn-2) = model.mantmparm.knots(2:end-1);    
    parprior(kk).vpvs = model.crustmparm.vpvs;    
    parprior(kk).cxi = model.crustmparm.xi;    
    parprior(kk).mxi = model.mantmparm.xi;    
end
Npass = sum(passed);

parprior(~passed) = [];
parprior = parprior';

% put into 1-D prior struct
prior.Npass = Npass;
prior.zsed = [parprior.zsed]';
prior.zmoh = [parprior.zmoh]';
prior.kcrust = [parprior.kcrust]';
prior.kmantle = [parprior.kmantle]';
prior.VSsedtop = [parprior.VSsedtop]';
prior.VSsedbot = [parprior.VSsedbot]';
prior.VScrusttop = [parprior.VScrusttop]';
prior.VScrustbot = [parprior.VScrustbot]';
prior.VSmantle = reshape([parprior.VSmantle],length(prior.zatdep),Npass)';
prior.fdVSsed = [parprior.fdVSsed]';
prior.fdVSmoh = [parprior.fdVSmoh]';
prior.vpvs = [parprior.vpvs]';
prior.cxi = [parprior.cxi]';
prior.mxi = [parprior.mxi]';
prior.zatdep = round(prior.zatdep);
prior.zmantle = prior.zatdep;

prior.Zkn_crust(1,:) = parprior(1).Zkn_crust;
prior.Zkn_mantle(1,:) = parprior(1).Zkn_mantle;
for ii = 2:Npass
prior.Zkn_crust(ii,:) = [parprior(ii).Zkn_crust,nan(1,par.mod.crust.kmax-length(parprior(ii).Zkn_crust))];
prior.Zkn_mantle(ii,:) = [parprior(ii).Zkn_mantle,nan(1,par.mod.mantle.kmax-length(parprior(ii).Zkn_mantle))];
end



% profile off
% profile viewer

end
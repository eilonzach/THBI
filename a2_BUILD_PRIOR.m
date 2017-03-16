function prior = a2_BUILD_PRIOR(par,Niter,zmantle)
% profile clear
% profile on

if nargin <2 || isempty(Niter)
    Niter = 1e5;
end

o = [];

prior = struct('Niter',Niter','Npass',0,...
               'zsed',o,'zmoh',o,...
               'kcrust',o,'kmantle',o,...
               'VSsedtop',o,'VSsedbot',o,'VScrusttop',o,'VScrustbot',o,...
               'VSmantle',repmat(o,1,6),'zmantle',zeros(6,1),...
               'fdVSsed',o,'fdVSmoh',o,'vpvs',o);

if nargin<3 || isempty(zmantle)
    prior.zmantle = linspace(par.mod.sed.hmax+par.mod.crust.hmax+0.1,par.mod.maxz,40)';
else
    prior.zmantle = zmantle;
end

% prep for par
parprior = prior;           
parprior(Niter) = prior(1); % set last val to give prior right dimension

passed = zeros(Niter,1);
           
for kk = 1:Niter
    model = b1_INITIATE_MODEL(par);
    ifpass = a1_TEST_CONDITIONS( model, par );
    
    if ~ifpass, continue; end
    
    passed(kk) = 1;
    
    parprior(kk).zsed = model.zsed; %#ok<PFOUS>
    parprior(kk).zmoh = model.zmoh;
    parprior(kk).VSsedtop = model.sedmparm.VS(1);
    parprior(kk).VSsedbot = model.sedmparm.VS(2);
    parprior(kk).VScrusttop = model.crustmparm.VS_sp(1);
    parprior(kk).VScrustbot = model.crustmparm.VS_sp(end);
    parprior(kk).VSmantle(1,:) = linterp(model.z,model.VS,prior.zmantle);
    parprior(kk).fdVSsed = model.fdVSsed;
    parprior(kk).fdVSmoh = model.fdVSmoh;
    parprior(kk).vpvs = model.crustmparm.vpvs;    
end
Npass = sum(passed);

parprior(~passed) = [];
parprior = parprior';

% put into 1-D prior struct
prior.Npass = Npass;
prior.zsed = [parprior.zsed]';
prior.zmoh = [parprior.zmoh]';
prior.VSsedtop = [parprior.VSsedtop]';
prior.VSsedbot = [parprior.VSsedbot]';
prior.VScrusttop = [parprior.VScrusttop]';
prior.VScrustbot = [parprior.VScrustbot]';
prior.VSmantle = reshape([parprior.VSmantle],length(prior.zmantle),Npass)';
prior.fdVSsed = [parprior.fdVSsed]';
prior.fdVSmoh = [parprior.fdVSmoh]';
prior.vpvs = [parprior.vpvs]';
prior.zmantle = round(prior.zmantle);



% profile off
% profile viewer

end
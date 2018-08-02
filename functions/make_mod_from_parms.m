function [ model ] = make_mod_from_parms( model,par )
% [ model ] = make_mod_from_parms( model )
% 
% make the 1D model vectors and important model values from the three
% structures sedparm, crustparm, and mantparm that actually house the
% important values that go towards building the model.

mps = model.sedmparm;
mpc = model.crustmparm;
mpm = model.mantmparm;

%% SEDIMENTS
zs = [0,mps.h]';

vs_sed = mps.VS(:);
vp_sed = sed_vs2vp(vs_sed);
rho_sed = sed_vs2rho(vs_sed);
xi_sed = ones(size(zs));

if diff(zs)==0; zs=[]; vs_sed=[]; vp_sed=[]; rho_sed=[]; xi_sed=[]; end

%% CRUST
cminz = mps.h;
cmaxz = mps.h+mpc.h;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';

vs_crust = sum(mpc.splines*diag(mpc.VS_sp),2);
vp_crust = vs_crust*mpc.vpvs;
rho_crust = sed_vs2rho(vs_crust); % use same expression as for sed (assume not ultramafic or other poorly-fit rocks)
xi_crust = mpc.xi*ones(size(zc));


%% MANTLE
mminz = cmaxz;
mmaxz = par.mod.maxz + model.selev; 
zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';

vs_mantle = sum(mpm.splines*diag(mpm.VS_sp),2);
vp_mantle = mantle_vs2vp(vs_mantle,zm );
rho_mantle = mantle_vs2rho(vs_mantle,zm );
xi_mantle = mpm.xi*ones(size(zm));


%% COLLATE
zz = [zs;zc;zm];
zz0 = zz-model.selev; % true depth, from sea level/ref. ellipsoid
Nz = length(zz);
zsed = mps.h;
zmoh = mps.h+mpc.h;
vs = [vs_sed;vs_crust;vs_mantle];
vp = [vp_sed;vp_crust;vp_mantle];
rho = [rho_sed;rho_crust;rho_mantle];
xi = [xi_sed;xi_crust;xi_mantle];

%% OUTPUT
model.z    = zz;
model.z0   = zz0;
model.VS   = vs;
model.VP   = vp;
model.rho  = rho;
model.Nz   = Nz;
model.zsed = zsed;
model.zmoh = zmoh;
model.vpvs = mpc.vpvs;
model.cxi  = mpc.xi;
model.mxi  = mpm.xi;
model.fdVSsed = 100*diff(vs(zz==zsed))./mean(vs(zz==zsed)); if zsed==0, model.fdVSsed = nan; end
model.fdVSmoh = 100*diff(vs(zz==zmoh))./mean(vs(zz==zmoh));
model.Sanis = 100*(xi-1); % in percentage 
model.Panis = zeros(Nz,1);

% re-order fields of model structure 
forder ={'z','z0','VS','VP','rho','Nz','zsed','zmoh','vpvs','cxi','mxi','selev',...
         'fdVSsed','fdVSmoh','Sanis','Panis',...
         'sedmparm','crustmparm','mantmparm','datahparm','M'};
     
model = orderfields(model,forder);

end


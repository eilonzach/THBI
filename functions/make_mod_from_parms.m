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

%% CRUST
cminz = mps.h;
cmaxz = mps.h+mpc.h;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';

vs_crust = sum(mpc.splines*diag(mpc.VS_sp),2);
vp_crust = vs_crust*mpc.vpvs;
rho_crust = sed_vs2rho(vs_crust); % use same expression as for sed (assume not ultramafic or other poorly-fit rocks)


%% MANTLE
mminz = cmaxz;
mmaxz = par.mod.maxz + model.selev; 
zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';

vs_mantle = sum(mpm.splines*diag(mpm.VS_sp),2);
vp_mantle = mantle_vs2vp(vs_mantle,zm );
rho_mantle = mantle_vs2rho(vs_mantle,zm );


%% COLLATE
zz = [zs;zc;zm];
zz0 = zz-model.selev; % true depth, from sea level/ref. ellipsoid
Nz = length(zz);
zsed = mps.h;
zmoh = mps.h+mpc.h;
vs = [vs_sed;vs_crust;vs_mantle];
vp = [vp_sed;vp_crust;vp_mantle];
rho = [rho_sed;rho_crust;rho_mantle];

%% OUTPUT
model.z    = zz;
model.z0   = zz0;
model.VS   = vs;
model.VP   = vp;
model.rho  = rho;
model.Nz   = Nz;
model.zsed = zsed;
model.zmoh = zmoh;
model.fdVSsed = 100*diff(vs(zz==zsed))./mean(vs(zz==zsed)); if zsed==0, model.fdVSsed = nan; end
model.fdVSmoh = 100*diff(vs(zz==zmoh))./mean(vs(zz==zmoh));
model.Sanis = zeros(Nz,1);
model.Panis = zeros(Nz,1);

% re-order fields of model structure 
forder ={'z','z0','VS','VP','rho','Nz','zsed','zmoh','selev',...
         'fdVSsed','fdVSmoh','Sanis','Panis',...
         'sedmparm','crustmparm','mantmparm','datahparm','M'};
     
model = orderfields(model,forder);

end


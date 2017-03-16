% model0 = model;
model = model0;

laymod = 1;

if laymod==1;
%% CRUST
dh = +5;
model.sedmparm.h = model.sedmparm.h + dh;
model.crustmparm.h = model.crustmparm.h - dh;

% modify splines in crust
cminz = model.sedmparm.h;
cmaxz = cminz + model.crustmparm.h;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
par.mod.dzsp = (cmaxz-cminz)/(model.crustmparm.Nkn-2);
cknots = [repmat(cminz,1,3),cminz:par.mod.dzsp:cmaxz,repmat(cmaxz,1,3)]';

iczt = find(model.z==model.zsed,1,'last');
iczb = find(model.z==model.zmoh,1,'first');
if dh>0
	zci = model.z(iczt:iczb);
    vci = model.VS(iczt:iczb);
elseif dh<0
	zci = [model.zsed+dh;model.z(iczt:iczb)];
    vci = [model.VS(iczt);model.VS(iczt:iczb);];
end

sp = fastBSpline.lsqspline(cknots,2,zci,vci); % interpolate onto current model
spbasis = sp.getBasis(zc); 
model.crustmparm.splines = spbasis(:,2:end-1);                
model.crustmparm.VS_kn = sp.weights(2:end-1); % pull back out spline coeff's from the interpolation        

elseif laymod==2

%% MANTLE
dh = -1;
model.crustmparm.h = model.crustmparm.h + dh;
par.mod.dz = 5;

% modify splines in crust
cminz = model.sedmparm.h;
cmaxz = cminz + model.crustmparm.h;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
par.mod.dzsp = (cmaxz-cminz)/(model.crustmparm.Nkn-2);
cknots = [repmat(cminz,1,3),cminz:par.mod.dzsp:cmaxz,repmat(cmaxz,1,3)]';

iczt = find(model.z==model.zsed,1,'last');
iczb = find(model.z==model.zmoh,1,'first');
if dh>0
	zci = [model.z(iczt:iczb);model.zmoh+dh];
    vci = [model.VS(iczt:iczb);model.VS(iczb)];
elseif dh<0
	zci = model.z(iczt:iczb);
    vci = model.VS(iczt:iczb);
end

sp = fastBSpline.lsqspline(cknots,2,zci,vci); % interpolate onto current model
spbasis = sp.getBasis(zc); 
model.crustmparm.splines = spbasis(:,2:end-1);                
model.crustmparm.VS_kn = sp.weights(2:end-1); % pull back out spline coeff's from the interpolation        

% modify splines in mantle
mminz = model.sedmparm.h + model.crustmparm.h;
mmaxz = par.mod.maxz + model.selev;
zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';
par.mod.dzsp = (mmaxz-mminz)/(model.mantmparm.Nkn-2);
mknots = [repmat(mminz,1,3),mminz:par.mod.dzsp:mmaxz,repmat(mmaxz,1,3)]';

imzt = find(model.z==model.zmoh,1,'last');
if dh>0
    zmi = model.z(imzt:end);
    vmi = model.VS(imzt:end);
elseif dh<0
     zmi = [model.zmoh+dh;model.z(imzt:end)];
     vmi = [model.VS(imzt);model.VS(imzt:end)];
end
sp = fastBSpline.lsqspline(mknots,2,zmi,vmi); % interpolate onto current model
spbasis = sp.getBasis(zm); 
model.mantmparm.splines = spbasis(:,2:end-1);
model.mantmparm.VS_kn = sp.weights(2:end-1); % pull back out spline coeff's from the interpolation


end

%% make etc
[ model ] = make_mod_from_parms( model,par );

plot_MOD_TRUEvsTRIAL(model0,model)
plot_PARAMETERISATION(model)

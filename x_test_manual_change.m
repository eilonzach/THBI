%% take current model, perturb it manually, see if fits improve
% profile clear
% profile on

%% current model data and misfit
% [ model ] = make_mod_from_parms( model,par );
cdata = b3_FORWARD_MODEL( model,Kbase,par,trudata,id,0);
cdata = predat_process(cdata,'PsRF',par);
cdata = predat_process(cdata,'SpRF',par);
cdata = predat_process(cdata,'PsRF_lo',par);
cdata = predat_process(cdata,'SpRF_lo',par);
% [cdata.SW.phV] = run_mineos(model,trudata.SW.periods,'cmodmake');
[ cmisfit ] = b4_CALC_MISFIT( trudata,cdata,par,0,SWwt ); % misfit has structures of summed errors
[ clog_likelihood,cmisfit ] = b5_CALC_LIKELIHOOD( cmisfit,trudata,model.datahparm,par );
    
%% perturb model
model2 = model;
% =================  enact perturbation  ================= 
% == change velocities
model2.mantmparm.VS_sp(1) = model2.mantmparm.VS_sp(1)+.0; 
model2.mantmparm.VS_sp(3) = model2.mantmparm.VS_sp(3)+.0; 
model2.mantmparm.VS_sp(4) = model2.mantmparm.VS_sp(4)-.0; 
model2.mantmparm.VS_sp(5) = model2.mantmparm.VS_sp(5)-.0; 
% model2.mantmparm.VS_sp(9) = TRUEmodel.mantmparm.VS_sp(9); 
model2.mantmparm.VS_sp(10) = model2.mantmparm.VS_sp(10)+.0; 
model2.crustmparm.VS_sp(1) = model2.crustmparm.VS_sp(1)-0.; 
% model2.crustmparm.VS_sp(2:4) = 3.7; 
% model2.mantmparm.VS_sp(:) = model2.mantmparm.VS_sp(1) ;
% model2.mantmparm.VS_sp(4) = model2.mantmparm.VS_sp(4)+0.2; 
% model2.mantmparm.VS_sp(6) = model2.mantmparm.VS_sp(6)-0.2; 

% == move moho
% model2.sedmparm.h = model2.sedmparm.h + 0.0000;
% model2.sedmparm.h = 2.9;
dhmoh = +.0001;
model2.crustmparm.h = model2.crustmparm.h + dhmoh;
% model2.crustmparm.h = 21.196;
if dhmoh~=0
% modify splines in crust
cminz = model2.sedmparm.h;
cmaxz = cminz + model2.crustmparm.h;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
par.mod.dzsp = (cmaxz-cminz)/(model2.crustmparm.Nkn-2);
cknots = [repmat(cminz,1,3),cminz:par.mod.dzsp:cmaxz,repmat(cmaxz,1,3)]';

iczt = find(model2.z==model2.zsed,1,'last');
iczb = find(model2.z==model2.zmoh,1,'first');
if dhmoh>0
	zci = [model2.z(iczt:iczb);model2.zmoh+dhmoh];
    vci = [model2.VS(iczt:iczb);model2.VS(iczb)];
elseif dhmoh<0
	zci = model2.z(iczt:iczb);
    vci = model2.VS(iczt:iczb);
end

sp = fastBSpline.lsqspline(cknots,2,zci,vci); % interpolate onto current model
spbasis = sp.getBasis(zc); 
model2.crustmparm.splines = spbasis(:,2:end-1);                
model2.crustmparm.VS_sp = sp.weights(2:end-1); % pull back out spline coeff's from the interpolation        

% modify splines in mantle
% mminz = model2.sedmparm.h + model2.crustmparm.h;
% mmaxz = par.mod.maxz + model2.selev;
% zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';
% par.mod.dzsp = (mmaxz-mminz)/(model2.mantmparm.Nkn-2);
% mknots = [repmat(mminz,1,3),mminz:par.mod.dzsp:mmaxz,repmat(mmaxz,1,3)]';

imzt = find(model2.z==model2.zmoh,1,'last');
if dhmoh>0
    zmi = model2.z(imzt:end);
    vmi = model2.VS(imzt:end);
elseif dhmoh<0
     zmi = [model2.zmoh+dhmoh;model2.z(imzt:end)];
     vmi = [model2.VS(imzt);model2.VS(imzt:end)];
end
% sp = fastBSpline.lsqspline(mknots,2,zmi,vmi); % interpolate onto current model
% spbasis = sp.getBasis(zm); 
% model2.mantmparm.splines = spbasis(:,2:end-1);
% model2.mantmparm.VS_sp = sp.weights(2:end-1); % pull back out spline coeff's from the interpolation
% 
fknots = linspace(0,1,9); 
% fknots(8) = 0.9*fknots(8);

[ model2.mantmparm.splines,model2.mantmparm.VS_sp,model2.mantmparm.knots ] = make_splines( model2,par,fknots,'mantle',zmi,vmi )

end
   
% == change Moho discontinuity
vt = model2.crustmparm.VS_sp(end);
vb = model2.mantmparm.VS_sp(1);
avV = (vt+vb)/2; % mean V of disc
dV0 = vb-vt; % V jump at disc
ddV = -0.;                  
davV = -0.0;                  
model2.crustmparm.VS_sp(end) = (avV+davV) - (dV0+ddV)/2;
model2.mantmparm.VS_sp(1)    = (avV+davV) + (dV0+ddV)/2;

% == change sed/crust discontinuity
vt = model2.sedmparm.VS(end);
vb = model2.crustmparm.VS_sp(1);
avV = (vt+vb)/2; % mean V of disc
dV0 = vb-vt; % V jump at disc
ddV = -0.0;                  
davV = 0;                  
model2.sedmparm.VS(end) = (avV+davV) - (dV0+ddV)/2;
model2.crustmparm.VS_sp(1)    = (avV+davV) + (dV0+ddV)/2;

% == change hyperparameters
model2.datahparm.sigmaPsRF = model2.datahparm.sigmaPsRF*1;

% =================  re-calc model  ================= 
[ model2 ] = make_mod_from_parms( model2,par );

% compute perturbation
modptb = calc_Vperturbation(model,model2);
ptbnorm = norm(modptb.dvsv) + norm(modptb.dvsh);

%% perturbed model data and misfit
predata = b3_FORWARD_MODEL( model2,Kbase,par,trudata,id,0);
predata = predat_process(predata,'PsRF',par);
predata = predat_process(predata,'SpRF',par);
predata = predat_process(predata,'PsRF_lo',par);
predata = predat_process(predata,'SpRF_lo',par);
% [predata.SW.phV] = run_mineos(model2,trudata.SW.periods,'cmodmake');
[ misfit ] = b4_CALC_MISFIT( trudata,predata,par,0,SWwt ); % misfit has structures of summed errors
[ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,trudata,model2.datahparm,par );
 



%% plot true, current, and trial models
plot_MOD_TRUEvsTRIAL( TRUEmodel, model )
subplot(1,2,1), hold on, h(1) = plot(model2.VS,model2.z,'b--');
subplot(1,2,2), hold on, h(2) = plot(model2.VP,model2.z,'b--');

%% plot true, current, and trial data
axs = plot_TRUvsPRE( trudata,cdata); 
axes(axs(5)), hold on,plot(predata.PsRF.tt,predata.PsRF.ZRT(:,1),'c','linewidth',2),
axes(axs(6)), hold on,plot(predata.PsRF.tt,predata.PsRF.ZRT(:,2),'c','linewidth',2),
axes(axs(7)), hold on,plot(predata.SpRF.tt,predata.SpRF.ZRT(:,1),'c','linewidth',2),
axes(axs(8)), hold on,plot(predata.SpRF.tt,predata.SpRF.ZRT(:,2),'c','linewidth',2),
axes(axs(9)), hold on,plot(predata.SW.periods,predata.SW.phV,'c.-','linewidth',3,'markersize',40);

pause(0.01)
figure(85)

%% report changes in misfit and likelihood
fprintf('=============================\n')
fprintf('              model   model1\n')
fprintf('-----------------------------\n')
fprintf('PsRFmisfit %8.3f  %8.3f   %8.3f\n',cmisfit.PsRF,misfit.PsRF,cmisfit.PsRF-misfit.PsRF)
fprintf('SpRFmisfit %8.3f  %8.3f   %8.3f\n',cmisfit.SpRF,misfit.SpRF,cmisfit.SpRF-misfit.SpRF)
fprintf('PsRFmisflo %8.3f  %8.3f   %8.3f\n',cmisfit.PsRF_lo,misfit.PsRF_lo,cmisfit.PsRF_lo-misfit.PsRF_lo)
fprintf('SpRFmisflo %8.3f  %8.3f   %8.3f\n',cmisfit.SpRF_lo,misfit.SpRF_lo,cmisfit.SpRF_lo-misfit.SpRF_lo)
fprintf('SWmisfit   %8.3f  %8.3f   %8.3f\n',cmisfit.SW,misfit.SW,cmisfit.SW-misfit.SW)
fprintf('Ps chi2    %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_ps,misfit.chi2_ps,cmisfit.chi2_ps-misfit.chi2_ps)
fprintf('Sp chi2    %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_sp,misfit.chi2_sp,cmisfit.chi2_sp-misfit.chi2_sp)
fprintf('Ps_lo chi2 %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_ps_lo,misfit.chi2_ps_lo,cmisfit.chi2_ps_lo-misfit.chi2_ps_lo)
fprintf('Sp_lo chi2 %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_sp_lo,misfit.chi2_sp_lo,cmisfit.chi2_sp_lo-misfit.chi2_sp_lo)
fprintf('SW chi2    %8.3f  %8.3f   %8.3f\n',cmisfit.chi2_SW,misfit.chi2_SW,cmisfit.chi2_SW-misfit.chi2_SW)
fprintf('Like''hood  %8.3f  %8.3f   %8.3f\n',clog_likelihood,log_likelihood,log_likelihood-clog_likelihood)
% 
% profile viewer
return
%% loop over hyperparamters
sigs = logspace(-3.5,0,60)';
l_hoods = zeros(length(sigs),3);
model1.datahparm = model.datahparm;
for iips = 1:length(sigs)
model1.datahparm.sigmaPsRF = sigs(iips);
[ l_hoods(iips,1) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaPsRF = model.datahparm.sigmaPsRF;
for iisp = 1:length(sigs)
model1.datahparm.sigmaSpRF = sigs(iisp);
[ l_hoods(iisp,2) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaSpRF = model.datahparm.sigmaSpRF;
for iisw = 1:length(sigs)
model1.datahparm.sigmaSW = sigs(iisw);
[ l_hoods(iisw,3) ] = b5_CALC_LIKELIHOOD( misfit,trudata,model1.datahparm,par );
end
model1.datahparm.sigmaSW = model.datahparm.sigmaSW;
figure(45),clf
loglog(sigs,l_hoods);
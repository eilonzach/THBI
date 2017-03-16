project = 'SYNTHETICS';


%% ------------------------- START ------------------------- 
run /Users/zeilon/Documents/MATLAB/BayesianJointInv/a0_STARTUP_BAYES;
cd([basedir,project]);
run parms/bayes_inv_parms

cp_ps = struct('samprate',par.synth.samprate,'pretime',-par.forc.Twin_Ps(1),...
               'prex',-par.forc.Twin_Ps(1),'postx',par.forc.Twin_Ps(2),'taperx',0.06,...
               'fhi',1/1,'flo',1/100,'npoles',2,'norm',0);
cp_sp = struct('samprate',par.synth.samprate,'pretime',-par.forc.Twin_Sp(1),...
               'prex',-par.forc.Twin_Sp(1),'postx',par.forc.Twin_Sp(2),'taperx',0.06,...
               'fhi',1/1,'flo',1/100,'npoles',2,'norm',0);


%% dz = 1 km
par.mod.dz=1;
tic
z0_SYNTH_MODEL_custommod(par,1); mod1 = TRUEmodel ; close(95)
toc
[ d1km ] = z1_SYNTH_DATA(par,0); % in ZRT format
toc
[ d1km.PsRF ] = predat_process( d1km.PsRF,cp_ps,par.mod.data.normdata,par.mod.data.decdata);
[ d1km.SpRF ] = predat_process( d1km.SpRF,cp_sp,par.mod.data.normdata,par.mod.data.decdata);

%% dz = 4 km
par.mod.dz=4;
tic
z0_SYNTH_MODEL_custommod(par,1); mod4 = TRUEmodel ; close(95)
toc
[ d4km ] = z1_SYNTH_DATA(par,0); % in ZRT format
toc
[ d4km.PsRF ] = predat_process( d4km.PsRF,cp_ps,par.mod.data.normdata,par.mod.data.decdata);
[ d4km.SpRF ] = predat_process( d4km.SpRF,cp_sp,par.mod.data.normdata,par.mod.data.decdata);

%% dz = 10 km
par.mod.dz=10;
tic
z0_SYNTH_MODEL_custommod(par,1); mod10 = TRUEmodel ; close(95)
toc
[ d10km ] = z1_SYNTH_DATA(par,0); % in ZRT format
toc
[ d10km.PsRF ] = predat_process( d10km.PsRF,cp_ps,par.mod.data.normdata,par.mod.data.decdata);
[ d10km.SpRF ] = predat_process( d10km.SpRF,cp_sp,par.mod.data.normdata,par.mod.data.decdata);



%% plot model
plot_MOD_TRUEvsTRIAL( mod1, mod4 )
subplot(1,2,1), hold on, h(1) = plot(mod10.VS,mod10.z,'b--');
subplot(1,2,2), hold on, h(2) = plot(mod10.VP,mod10.z,'b--');


%% plot data
plot_TRUvsPRE(d1km,d4km)
subplot(231), hold on,plot(d10km.PsRF.tt,d10km.PsRF.ZRT(:,1),'c','linewidth',2), ylabel('Z pre','fontsize',18)
subplot(234), hold on,plot(d10km.PsRF.tt,d10km.PsRF.ZRT(:,2),'c','linewidth',2), ylabel('R pre','fontsize',18)
subplot(232), hold on,plot(d10km.SpRF.tt,d10km.SpRF.ZRT(:,1),'c','linewidth',2),  ylabel('Z pre','fontsize',18)
subplot(235), hold on,plot(d10km.SpRF.tt,d10km.SpRF.ZRT(:,2),'c','linewidth',2), ylabel('R pre','fontsize',18)
subplot(2,3,[3,6]), hold on,plot(d10km.SW.periods,d10km.SW.phV,'c.-','linewidth',3,'markersize',40);

subplot(2,3,[3,6]), hl = legend('1km','4km','10km','location','southeast'); set(hl,'fontsize',15)

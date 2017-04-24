function [ data ] = z1_SYNTH_DATA(par,ifplot)

global TLM TRUEmodel 

% synthetic data parameters
samprate = par.synth.samprate;
SWperiods = par.synth.SWperiods;
gcarc = 60;
tt = tauptime('deg',gcarc,'phases','P,S');
P_inc = rayp2inc(tt(1).rayparameter,TLM.Vp(end),6371-TLM.zlayb(end));
S_inc = rayp2inc(tt(2).rayparameter,TLM.Vs(end),6371-TLM.zlayb(end));
%% ===================  CALCULATE RECEIVER FUNCTIONS  ===================

[trudat_ps,tt_ps]  = run_propmat(TLM,'TegPs','Ps',samprate, P_inc, par.forc.synthperiod); % in R,T,Z
[trudat_sp,tt_sp] = run_propmat(TLM,'TegSp','Sp',samprate, S_inc, par.forc.synthperiod); % in R,T,Z
% channel order
trudat_ps = trudat_ps(:,[3,1,2]); % in Z,R,T
trudat_sp = trudat_sp(:,[3,1,2]); % in Z,R,T
    
tt_Par = mean(tt_ps(trudat_ps(:,1)==min(trudat_ps(:,1))));% estimate main P-wave arrival time from first big downswing
tt_Sar = mean(tt_sp(trudat_sp(:,2)==min(trudat_sp(:,2)))); % estimate main S-wave arrival time from first big downswing

tt_ps = tt_ps - tt_Par;
tt_sp = tt_sp - tt_Sar;
tt_ps = round_level(tt_ps,0.001);
tt_sp = round_level(tt_sp,0.001);

Ps_widewind = [-10 50];
Sp_widewind = [-50 10];


inwind_ps = (tt_ps >= Ps_widewind(1)) & (tt_ps < Ps_widewind(2)); 
inwind_sp = (tt_sp >= Sp_widewind(1)) & (tt_sp < Sp_widewind(2)); 

% crop
trudat_ps = trudat_ps(inwind_ps,:);
trudat_sp = trudat_sp(inwind_sp,:);
tt_ps = tt_ps(inwind_ps);
tt_sp = tt_sp(inwind_sp);

% taper
trudat_ps = flat_hanning_win(tt_ps,trudat_ps,Ps_widewind(1),Ps_widewind(2),3); % 3s taper
trudat_sp = flat_hanning_win(tt_sp,trudat_sp,Sp_widewind(1),Sp_widewind(2),3); % 3s taper

% normalise to small energy
normf_ps = trudat_ps(:,1)'*trudat_ps(:,1) + trudat_ps(:,2)'*trudat_ps(:,2);
trudat_ps = trudat_ps/sqrt(normf_ps/2);

normf_sp = trudat_sp(:,1)'*trudat_sp(:,1) + trudat_sp(:,2)'*trudat_sp(:,2);
trudat_sp = trudat_sp/sqrt(normf_sp/2);

% add noise
trudat_ps = trudat_ps + random('norm',0,par.synth.noise_sigmaPsRF,size(trudat_ps));
trudat_sp = trudat_sp + random('norm',0,par.synth.noise_sigmaSpRF,size(trudat_sp)); 

if ifplot
	figure(58),clf,set(gcf,'pos',[2 275 1047 830])
    %Ps
    subplot(411),plot(tt_ps,trudat_ps(:,1),'k','linewidth',2), xlim(par.datprocess.Twin.PsRF), ylabel('Z tru','fontsize',18), title('Ps TRUE','fontsize',22)
    subplot(412),plot(tt_ps,trudat_ps(:,2),'r','linewidth',2), xlim(par.datprocess.Twin.PsRF), ylabel('R tru','fontsize',18)
    %Sp
    subplot(413),plot(tt_sp,trudat_sp(:,1),'k','linewidth',2), xlim(par.datprocess.Twin.SpRF), ylabel('Z tru','fontsize',18), title('Sp TRUE','fontsize',22)
    subplot(414),plot(tt_sp,trudat_sp(:,2),'r','linewidth',2), xlim(par.datprocess.Twin.SpRF), ylabel('R tru','fontsize',18)
end

PsRF = struct('ZRT',trudat_ps,'tt',tt_ps,'rayp',tt(1).rayparameter,'inc',P_inc,'samprate',samprate,'nsamp',length(tt_ps),'sigma',par.mod.data.prior_sigmaPsRF);
SpRF = struct('ZRT',trudat_sp,'tt',tt_sp,'rayp',tt(2).rayparameter,'inc',S_inc,'samprate',samprate,'nsamp',length(tt_sp),'sigma',par.mod.data.prior_sigmaSpRF);



%% ===================  CALCULATE PHASE VELOCITIES  ===================
% % Use Menke phV solver method on layered model 
% phinmod = [TLM.zlayb-TLM.zlayt,TLM.Vp,TLM.Vs,TLM.rho];
% trudat_cph = Calc_Ray_dispersion(SWperiods,phinmod,1,2000,0);
% % Use MINEOS
[trudat_cph] = run_mineos(TRUEmodel,SWperiods,0,'initmod',0);

% add noise
trudat_cph = trudat_cph + random('norm',0,par.synth.noise_sigmaSW,size(trudat_cph));

if ifplot
    figure(13);plot(SWperiods,trudat_cph,'o')
end

SW = struct('periods',SWperiods,'phV',trudat_cph,'sigma',par.mod.data.prior_sigmaSW);

%% ===================  MAKE DATA STRUCTURE  ===================
data = struct('PsRF',PsRF,'SpRF',SpRF,'SW',SW);




end


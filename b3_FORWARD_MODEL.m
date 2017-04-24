function [ predata ] = b3_FORWARD_MODEL( model,Kbase,par,data,ID,ifplot )
%[ predat_ps,predat_sp ] = b3_FORWARD_MODEL( model,par,ID,ifplot )
% 
%   Do forward model to calculate predicted data.
% 
% INPUTS
%   model   - model structure
%   Kbase   - structure with kernel model, depth kernels, and its phase vels
%   par     - parameters structure
%   data    - obs data structure with all datatypes 
%   ID      - unique ID for the propmat script to avoid overwriting files
%             if running in parallel.
%   ifplot  - flag with option to plot (1) or not (0)
% 
% OUTPUTS
%   predata - structure identical to input data structure, but with
%             predicted data, rather than observed data
%%
% An important component is the layerising of the model - conversion of
% continuous model into a bunch of layers, with coarseness partly
% determined by the minimum dVs in any layer (specified as an input). The
% layerised 1D model is also output from this function.

%% ===================  PREPARE DATA STRUCTURE  ===================
predata = data;



%% ===================  LAYERISE PROFILE  ===================
[zlayt,zlayb,Vslay,Vplay,rholay] = ...
    layerise(model.z,model.VS,par.forc.mindV,0,model.VP,model.rho); 
nlay = length(Vslay);
laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay);
if any(isnan(laymodel.rho))
    error('NaN densities')
end

if ifplot
    figure(1); clf, hold on
    plot(model.VS,model.z,'-ko')
    plot(model.VP,model.z,'-ko')
    zlayp = reshape([laymodel.zlayt';laymodel.zlayb'],2*laymodel.nlay,1);
    vslayp = reshape([laymodel.Vs';laymodel.Vs'],2*laymodel.nlay,1);
    vplayp = reshape([laymodel.Vp';laymodel.Vp'],2*laymodel.nlay,1);
    plot(vslayp,zlayp,'-ro')
    plot(vplayp,zlayp,'-ro')
    set(gca,'ydir','reverse','ylim',[0, max(model.z)],'xlim',[0.9*min(model.VS) 1.1*max(model.VP)])
    set(gcf,'pos',[41   282   729   823]);
end

%% ===================  PS RFs FROM PROPAGATOR MATRIX  ====================

if any(strcmp(par.inv.datatypes,'PsRF'))
    [ unique_rayps_P,irayps_P ] = rayp_vals( [data.PsRF.rayp] );
    for ir = 1:length(unique_rayps_P)
        rayp = unique_rayps_P(ir);
        samprate = unique([data.PsRF(irayps_P==ir).samprate]);
        P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
        [predat_ps,tt_ps] = run_propmat(laymodel,ID,'Ps',samprate, P_inc, par.forc.synthperiod,par.forc.nsamps);
        predat_ps = predat_ps(:,[3,1,2]); % in Z,R,T
        tt_ps_Par = mean(tt_ps(predat_ps(:,1)==min(predat_ps(:,1))));% estimate main P-wave arrival time from first big downswing
        tt_ps = tt_ps - tt_ps_Par;
        tt_ps = round_level(tt_ps,0.001);

        Ps_widewind = [-10 50];
        inwind = (tt_ps >= Ps_widewind(1)) & (tt_ps < Ps_widewind(2)); 
        % crop
        predat_ps = predat_ps(inwind,:);
        tt_ps = tt_ps(inwind);
        if ~any(inwind)
            tt_ps = []; predat_ps = nan; % NO GOOD DATA
        else
            % taper
            predat_ps = flat_hanning_win(tt_ps,predat_ps,Ps_widewind(1),Ps_widewind(2),1); % 1s taper
            % normalise to unit energy
            % normf_ps = predat_ps(:,1)'*predat_ps(:,1) + predat_ps(:,2)'*predat_ps(:,2) + predat_ps(:,3)'*predat_ps(:,3);
            % predat_ps = predat_ps/sqrt(normf_ps);'
        end
        
        % -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_P==ir);
        for iir = 1:length(inds)
            predata.PsRF(inds(iir),1).ZRT=predat_ps;
            predata.PsRF(inds(iir),1).tt=tt_ps;
            predata.PsRF(inds(iir),1).nsamp = length(predata.PsRF(inds(iir)).ZRT);
        end
    end
else
    predata.PsRF = [];    
end

% Ps_lo
if any(strcmp(par.inv.datatypes,'PsRF_lo'))
    predata.PsRF_lo = predata.PsRF;
else 
    predata.PsRF_lo = [];
end


    
%% ===================  SP RFs FROM PROPAGATOR MATRIX  ====================
if any(strcmp(par.inv.datatypes,'SpRF'))
    [ unique_rayps_S,irayps_S ] = rayp_vals( [data.SpRF.rayp] );
    for ir = 1:length(unique_rayps_S)
        rayp = unique_rayps_S(ir);
        samprate = unique([data.SpRF(irayps_S==ir).samprate]);
        S_inc = rayp2inc(rayp,laymodel.Vs(end),6371-laymodel.zlayb(end));
        % check if inhomogeneous
        if isreal(asind(laymodel.Vp*sind(S_inc)./laymodel.Vs(end))) % 
            [predat_sp,tt_sp] = run_propmat(laymodel,ID,'Sp',samprate, S_inc, par.forc.synthperiod,par.forc.nsamps);
            predat_sp = predat_sp(:,[3,1,2]); % in Z,R,T
            tt_sp_Sar = mean(tt_sp(predat_sp(:,2)==min(predat_sp(:,2)))); % estimate main S-wave arrival time from first big downswing
            tt_sp = tt_sp - tt_sp_Sar;
            tt_sp = round_level(tt_sp,0.001);

            Sp_widewind = [-50 10];
            inwind = (tt_sp >= Sp_widewind(1)) & (tt_sp < Sp_widewind(2)); 
            % crop
            predat_sp = predat_sp(inwind,:);
            tt_sp = tt_sp(inwind);

            if ~any(inwind)
                tt_ps = []; predat_ps = nan; % NO GOOD DATA
            else
                % taper
                predat_sp = flat_hanning_win(tt_sp,predat_sp,Sp_widewind(1),Sp_widewind(2),1); % 1s taper
                % normalise to unit energy
                % normf_sp = predat_sp(:,1)'*predat_sp(:,1) + predat_sp(:,2)'*predat_sp(:,2) + predat_sp(:,3)'*predat_sp(:,3);
                % predat_sp = predat_sp/sqrt(normf_sp);
            end
        else
            tt_sp = NaN; predat_sp = NaN;    
        end
        
% -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_S==ir);
        for iir = 1:length(inds)
            predata.SpRF(inds(iir),1).ZRT=predat_sp;
            predata.SpRF(inds(iir),1).tt=tt_sp;
            predata.SpRF(inds(iir),1).nsamp = length(predata.SpRF(inds(iir)).ZRT);
        end
    end
else
    predata.SpRF = [];
end

% Sp_lo
if any(strcmp(par.inv.datatypes,'SpRF_lo'))
    predata.SpRF_lo = predata.SpRF;
else 
    predata.SpRF_lo = [];
end


%% ifplot....

if ifplot
    figure(2); clf, hold on
    comps = {'VERTICAL','RADIAL','TRANSVERSE'}; 
    for ip = 1:3
    subplot(3,2,2*ip-1)
    plot(tt_ps,predat_ps(:,ip),'Linewidth',2)
    xlim(par.datprocess.Twin_Ps);
    ylabel(comps{ip},'fontsize',19,'fontweight','bold')
    if ip == 1, title('Ps','fontsize',22,'fontweight','bold'), end
    subplot(3,2,2*ip)
    plot(tt_sp,predat_sp(:,ip),'Linewidth',2)
    xlim(par.datprocess.Twin.SpRF);
    if ip == 1, title('Sp','fontsize',22,'fontweight','bold'), end
    end
    set(gcf,'position',[780         282        1058         823]);
end % on ifplot

%% ===================  CALCULATE PHASE VELOCITIES  ===================
if any(strcmp(par.inv.datatypes,'SW'))
% phinmod = [laymodel.zlayb-laymodel.zlayt,laymodel.Vp,laymodel.Vs,laymodel.rho];
% predat_cph = Calc_Ray_dispersion(data.SW.periods,phinmod,1,2000,0);

% calc. perturbation values from 0==>1 and use to calc dc/c
[ modptb ] = calc_Vperturbation(Kbase.modelk,model);
Np = length(Kbase.phV);
est_dc_c = zeros(Np,1);
kfld = {'Vsv','Vsh','Vpv','Vph','rho';'dvsv','dvsh','dvpv','dvph','drho'};
for ip = 1:Np
    for ik = 1:size(kfld,2)
        dr = diff(modptb.Z).*1e3;
        K = midpts(Kbase.K{ip}.(kfld{1,ik}));
        dval = midpts(modptb.(kfld{2,ik}));
        est_dc_c(ip) = est_dc_c(ip) + dr'*(K.*dval);
    end
end
% estimated model1 phase velocities
predat_cph = (1+est_dc_c).*Kbase.phV;


% figure(13);plot(surfperiods,predat_cph,'o')
else
   predat_cph = [];
end

% -----------------  PUT INTO DATA STRUCTURE  -----------------
% SW
if any(strcmp(par.inv.datatypes,'SW'))
    predata.SW.phV = predat_cph;
end


% %% ===================  PUT INTO DATA STRUCTURE  ===================
% predata = data;
% 
% % Ps
% if any(strcmp(par.inv.datatypes,'PsRF'))
% predata.PsRF = predata.PsRF(1);
% predata.PsRF.ZRT = predat_ps;
% predata.PsRF.tt = tt_ps;
% predata.PsRF.samprate = data.PsRF.samprate;
% predata.PsRF.nsamp = length(predat_ps);
% end
% 
% % Sp
% if any(strcmp(par.inv.datatypes,'SpRF'))
% predata.SpRF = predata.SpRF(1);
% predata.SpRF.ZRT = predat_sp;
% predata.SpRF.tt = tt_sp;
% predata.SpRF.samprate = data.SpRF.samprate;
% predata.SpRF.nsamp = length(predat_sp);
% end
% 
% % SW
% if any(strcmp(par.inv.datatypes,'SW'))
% predata.SW.phV = predat_cph;
% end
% 
% % Ps_lo
% if any(strcmp(par.inv.datatypes,'PsRF_lo'))
% predata.PsRF_lo = predata.PsRF_lo(1);
% predata.PsRF_lo.ZRT = predat_ps;
% predata.PsRF_lo.tt = tt_ps;
% predata.PsRF_lo.samprate = data.PsRF.samprate;
% predata.PsRF_lo.nsamp = length(predat_ps);
% end
% 
% % Sp_lo
% if any(strcmp(par.inv.datatypes,'SpRF_lo'))
% predata.SpRF_lo = predata.SpRF_lo(1);
% predata.SpRF_lo.ZRT = predat_sp;
% predata.SpRF_lo.tt = tt_sp;
% predata.SpRF_lo.samprate = data.SpRF.samprate;
% predata.SpRF_lo.nsamp = length(predat_sp);
% end

end


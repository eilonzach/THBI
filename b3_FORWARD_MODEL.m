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

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

% [predata.PsRF.Vp_surf] = deal(mean([predata.PsRF.Vp_surf]));
% [predata.PsRF.Vs_surf] = deal(mean([predata.PsRF.Vs_surf]));


%% ===================  LAYERISE PROFILE  ===================
[zlayt,zlayb,Vslay] = ...
    layerise(model.z,model.VS,par.forc.mindV,0); 
nlay = length(Vslay);

% S to P and rho structure
xs = 1:find(zlayb==model.zsed); if model.zsed ==0, xs = []; end
xc = find(zlayt==model.zsed):find(zlayb==model.zmoh);
xm = find(zlayt==model.zmoh):nlay;
Vplay = [sed_vs2vp(Vslay(xs));...
         model.crustmparm.vpvs*Vslay(xc);...
         mantle_vs2vp(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
rholay = [sed_vs2rho(Vslay([xs,xc]));...
          mantle_vs2rho(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
xilay = [zeros(length(xs),1);...
         model.crustmparm.xi*ones(length(xc),1);...
         model.mantmparm.xi*ones(length(xm),1)]; % S radial anisotropy
philay = ones(nlay,1); % P radial anisotropy
etalay = ones(nlay,1); % eta anisotropy

laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay,'xi',xilay,'phi',philay,'eta',etalay);
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

if any(strcmp(pdtyps(:,2),'Ps'))
    Psdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Ps'),1,'first')};
    [ unique_rayps_P,irayps_P ] = rayp_vals( [predata.(Psdat).rayp] );
    for ir = 1:length(unique_rayps_P)
        rayp = unique_rayps_P(ir);
        samprate = unique([predata.(Psdat)(irayps_P==ir).samprate]);
        P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
        [predat_ps,tt_ps] = run_propmat(laymodel,ID,'Ps',samprate, P_inc, par.forc.synthperiod,par.forc.nsamps);
        % pad with zeros
        tt_ps = [tt_ps(1) + [-1000:-1]'./samprate; tt_ps ;tt_ps(end) + [1:1000]'./samprate];
        predat_ps = [zeros(1000,3);predat_ps;zeros(1000,3)];
        %correct corrdinate order
        predat_ps_ZRT = predat_ps(:,[3,1,2]); % in Z,R,T
        if strcmp(par.forc.PSVorZR,'PSV')
            clear predat_ps_PSV;
            % convert to P, SV
            [predat_ps_PSV(:,1),predat_ps_PSV(:,2)] = ...
                Rotate_XZ_to_PSV(predat_ps_ZRT(:,2),-predat_ps_ZRT(:,1),...
                mean([predata.(Psdat).Vp_surf]),mean([predata.(Psdat).Vs_surf]),...
                rayp_sdeg2skm(rayp,laymodel.zlayb(end)));
        elseif strcmp(par.forc.PSVorZR,'ZR')
            % keep as ZR (but kill T; Z positive UP)
            predat_ps_PSV = predat_ps_ZRT(:,[1,2]);         
        end
        predat_ps_PSV = predat_ps_PSV./maxab(predat_ps_PSV(:,1)); % normalise on parental max, make positive
        tt_ps_Par = mean(tt_ps(predat_ps_PSV(:,1)==max(predat_ps_PSV(:,1))));% estimate main P-wave arrival time from first big upswing
        tt_ps = tt_ps - tt_ps_Par;
        tt_ps = round_level(tt_ps,0.001);

        Ps_widewind = [-10 50];
        inwind = (tt_ps >= Ps_widewind(1)) & (tt_ps < Ps_widewind(2)); 
        % crop
        predat_ps_PSV = predat_ps_PSV(inwind,:);
        tt_ps = tt_ps(inwind);
        if ~any(inwind)
            tt_ps = []; predat_ps_PSV = nan; % NO GOOD DATA
        else
            % taper
            predat_ps_PSV = flat_hanning_win(tt_ps,predat_ps_PSV,Ps_widewind(1),Ps_widewind(2),1); % 1s taper
            % normalise to unit energy
            % normf_ps = predat_ps(:,1)'*predat_ps(:,1) + predat_ps(:,2)'*predat_ps(:,2) + predat_ps(:,3)'*predat_ps(:,3);
            % predat_ps = predat_ps/sqrt(normf_ps);'
        end
        
        % -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_P==ir);
        for iir = 1:length(inds)
            predata.(Psdat)(inds(iir),1).PSV=predat_ps_PSV; 
            predata.(Psdat)(inds(iir),1).tt=tt_ps;
            predata.(Psdat)(inds(iir),1).nsamp = length(predata.(Psdat)(inds(iir)).PSV);
        end
    end
end

%% ===================  SP RFs FROM PROPAGATOR MATRIX  ====================
if any(strcmp(pdtyps(:,2),'Sp'))
    Spdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Sp'),1,'first')};
    [ unique_rayps_S,irayps_S ] = rayp_vals( [predata.(Spdat).rayp] );
    for ir = 1:length(unique_rayps_S)
        rayp = unique_rayps_S(ir);
        samprate = unique([predata.(Spdat)(irayps_S==ir).samprate]);
        S_inc = rayp2inc(rayp,laymodel.Vs(end),6371-laymodel.zlayb(end));
        P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
        % check if S inhomogeneous
        if isreal(asind(laymodel.Vs*sind(S_inc)./laymodel.Vs(end))) % 
            
            % find layers where S to P conversion will not go inhomogeneous
            Play_incs = asind(laymodel.Vp*sind(P_inc)./laymodel.Vp(end));
            if any(~isreal(Play_incs))
                nimagplay = [1:find(imag(Play_incs),1,'first')-1];
                fns = fieldnames(laymodel);
                laymodel_Suse = laymodel; 
                laymodel_Suse.nlay = length(nimagplay);
                for jj = 1:length(fns)
                    if length(laymodel.(fns{jj}))==1, continue; end
                    laymodel_Suse.(fns{jj}) = laymodel_Suse.(fns{jj})(nimagplay);
                end
                % set appropriate S_inc for the actual base
                S_inc = rayp2inc(rayp,laymodel_Suse.Vs(end),6371-laymodel_Suse.zlayb(end));
            else
                laymodel_Suse = laymodel;
            end
            [predat_sp,tt_sp] = run_propmat(laymodel_Suse,ID,'Sp',samprate, S_inc, par.forc.synthperiod,par.forc.nsamps);
            % pad with zeros
            tt_sp = [tt_sp(1) + [-1000:-1]'./samprate; tt_sp ;tt_sp(end) + [1:1000]'./samprate];
            predat_sp = [zeros(1000,3);predat_sp;zeros(1000,3)];
            %correct corrdinate order
            predat_sp_ZRT = predat_sp(:,[3,1,2]); % in Z,R,T
            if strcmp(par.forc.PSVorZR,'PSV')
                clear predat_sp_PSV;
                % convert to P, SV
                [predat_sp_PSV(:,1),predat_sp_PSV(:,2)] = ...
                    Rotate_XZ_to_PSV(predat_sp_ZRT(:,2),-predat_sp_ZRT(:,1),...
                    mean([predata.(Spdat).Vp_surf]),mean([predata.(Spdat).Vs_surf]),...
                    rayp_sdeg2skm(rayp,laymodel_Suse.zlayb(end)));
            elseif strcmp(par.forc.PSVorZR,'ZR')
                % keep as ZR (but kill T; Z positive UP)
                predat_sp_PSV = predat_sp_ZRT(:,[1,2]);         
            end
            predat_sp_PSV = predat_sp_PSV./maxab(predat_sp_PSV(:,2)); % normalise on parental max, make positive
            tt_sp_Sar = mean(tt_sp(predat_sp_PSV(:,2)==max(predat_sp_PSV(:,2)))); % estimate main S-wave arrival time from first big upswing
            tt_sp = tt_sp - tt_sp_Sar;
            tt_sp = round_level(tt_sp,0.001);

            Sp_widewind = [-50 10];
            inwind = (tt_sp >= Sp_widewind(1)) & (tt_sp < Sp_widewind(2)); 
            % crop
            predat_sp_PSV = predat_sp_PSV(inwind,:);
            tt_sp = tt_sp(inwind);

            if ~any(inwind)
                tt_ps = []; predat_sp_PSV = nan; % NO GOOD DATA
            else
                % taper
                predat_sp_PSV = flat_hanning_win(tt_sp,predat_sp_PSV,Sp_widewind(1),Sp_widewind(2),1); % 1s taper
                % normalise to unit energy
                % normf_sp = predat_sp(:,1)'*predat_sp(:,1) + predat_sp(:,2)'*predat_sp(:,2) + predat_sp(:,3)'*predat_sp(:,3);
                % predat_sp = predat_sp/sqrt(normf_sp);
            end
        else
            tt_sp = NaN; predat_sp_PSV = NaN;    
        end
        
% -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_S==ir);
        for iir = 1:length(inds)
            predata.(Spdat)(inds(iir),1).PSV=predat_sp_PSV; 
            predata.(Spdat)(inds(iir),1).tt=tt_sp;
            predata.(Spdat)(inds(iir),1).nsamp = length(predata.(Spdat)(inds(iir)).PSV);
        end
    end
end



%% distribute data for different processing (e.g. _lo, _cms)
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt}; 
    pdt = parse_dtype( dtype ); 
    if strcmp(pdt{1},'BW') && (~strcmp(pdt{3},'def') || ~strcmp(pdt{4},'def'))
        if any(strcmp(par.inv.datatypes,['BW_',pdt{2}])) % only if there IS a standard!
            predata.(dtype) = predata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
        end
    end
end


%% ifplot....

if ifplot
    figure(2); clf, hold on
    comps = {'VERTICAL','RADIAL','TRANSVERSE'}; 
    for ip = 1:2
    subplot(3,2,2*ip-1)
    plot(tt_ps,predat_ps_PSV(:,ip),'Linewidth',2)
    xlim(par.datprocess.Ps.Twin.def);
    ylabel(comps{ip},'fontsize',19,'fontweight','bold')
    if ip == 1, title('Ps','fontsize',22,'fontweight','bold'), end
    subplot(3,2,2*ip)
    plot(tt_sp,predat_sp_PSV(:,ip),'Linewidth',2)
    xlim(par.datprocess.Sp.Twin.def);
    if ip == 1, title('Sp','fontsize',22,'fontweight','bold'), end
    end
    set(gcf,'position',[780         282        1058         823]);
end % on ifplot

%% ===================  CALCULATE PHASE VELOCITIES  ===================
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    pdtyp = parse_dtype(dtype);
    if ~strcmp(pdtyp{1},'SW'), continue, end

    % calc. perturbation values from 0==>1 and use to calc dc/c
    [ modptb ] = calc_Vperturbation(Kbase.modelk,model);
    
    Np = length(Kbase.(pdtyp{2}).phV);
    est_dc_c = zeros(Np,1);
    
    kfld = {'Vsv','Vsh','Vpv','Vph','rho';'dvsv','dvsh','dvpv','dvph','drho'};
    
    zind = find(modptb.Z<max(Kbase.(pdtyp{2}).Kph{1}.Z/1000)); 

    for ip = 1:Np
        for ik = 1:size(kfld,2)
            dr = diff(modptb.Z(zind))*1e3; % km to m
            K = midpts(Kbase.(pdtyp{2}).(['K',pdtyp{3}(1:2)]){ip}.(kfld{1,ik})(zind));
            dval = midpts(modptb.(kfld{2,ik})(zind));
            est_dc_c(ip) = est_dc_c(ip) + dr'*(K.*dval);
        end
    end

    % estimated model1 phase velocities
    predata.(dtype).phV = (1+est_dc_c).*Kbase.(pdtyp{2}).(pdtyp{3});

end


end


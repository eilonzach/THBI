function [ final_predata ] = c5_FINAL_FORWARD_MODEL( final_model,par,data )
%[ final_predata ] = c5_FINAL_FORWARD_MODEL( final_model,par,data )

ID = 'finalmod';

%% ===================  PUT INTO DATA STRUCTURE  ===================
final_predata = data;

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

%% ===================  LAYERISE PROFILE  ===================
% [zlayt,zlayb,Vslay,Vplay,rholay] = ...
%     layerise(final_model.Z,final_model.VSav,par.forc.mindV,1,final_model.VPav,final_model.rhoav); 
% nlay = length(Vslay);
% laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay);
% if any(isnan(laymodel.rho))
%     error('NaN densities')
% end

[zlayt,zlayb,Vslay] = ...
    layerise(final_model.Z,final_model.VSav,par.forc.mindV,1); 
nlay = length(Vslay);

% S to P and rho structure
zsed = final_model.Zd(1).mu;
zmoh = final_model.Zd(2).mu;
xs = 1:mindex(zlayb,zsed); if zsed==0, xs = []; end
xc = mindex(zlayt,zsed):mindex(zlayb,zmoh);
xm = mindex(zlayt,zmoh):nlay;
Vplay = [sed_vs2vp(Vslay(xs));...
         final_model.vpvsav*Vslay(xc);...
         mantle_vs2vp(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
rholay = [sed_vs2rho(Vslay([xs,xc]));...
          mantle_vs2rho(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
xilay = [zeros(length(xs),1);...
         final_model.xicrav*ones(length(xc),1);...
         final_model.ximaav*ones(length(xm),1)]; % S radial anisotropy
philay = ones(nlay,1); % P radial anisotropy
etalay = ones(nlay,1); % eta anisotropy

laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay,'xi',xilay,'phi',philay,'eta',etalay);
if any(isnan(laymodel.rho))
    error('NaN densities')
end


%% ===================  PS RFs FROM PROPAGATOR MATRIX  ====================

if any(strcmp(pdtyps(:,2),'Ps'))
    Psdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Ps'),1,'first')};
    [ unique_rayps_P,irayps_P ] = rayp_vals( [final_predata.(Psdat).rayp] );
    for ir = 1:length(unique_rayps_P)
        rayp = unique_rayps_P(ir);
        samprate = unique([final_predata.(Psdat)(irayps_P==ir).samprate]);
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
                mean([final_predata.(Psdat).Vp_surf]),mean([final_predata.(Psdat).Vs_surf]),...
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
        % taper
        predat_ps_PSV = flat_hanning_win(tt_ps,predat_ps_PSV,Ps_widewind(1),Ps_widewind(2),3); % 3s taper
        % normalise to unit energy
        %normf_ps = predat_ps_PSV(:,1)'*predat_ps_PSV(:,1) + predat_ps_PSV(:,2)'*predat_ps_PSV(:,2);
        %predat_ps_PSV = predat_ps_PSV/sqrt(normf_ps);
        % -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_P==ir);
        for iir = 1:length(inds)
            final_predata.(Psdat)(inds(iir),1).PSV=predat_ps_PSV; 
            final_predata.(Psdat)(inds(iir),1).tt=tt_ps;
            final_predata.(Psdat)(inds(iir),1).nsamp = length(final_predata.(Psdat)(inds(iir)).PSV);
        end
    end
end

%% ===================  SP RFs FROM PROPAGATOR MATRIX  ====================
if any(strcmp(pdtyps(:,2),'Sp'))
if any(~strcmp(pdtyps(strcmp(pdtyps(:,2),'Sp'),3),'ccp'))
    Spdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Sp'),1,'first')};
    [ unique_rayps_S,irayps_S ] = rayp_vals( [final_predata.(Spdat).rayp] );
    for ir = 1:length(unique_rayps_S)
        rayp = unique_rayps_S(ir);
        samprate = unique([final_predata.(Spdat)(irayps_S==ir).samprate]);
        S_inc = rayp2inc(rayp,laymodel.Vs(end),6371-laymodel.zlayb(end));
        
        P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
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
                mean([final_predata.(Spdat).Vp_surf]),mean([final_predata.(Spdat).Vs_surf]),...
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
        % taper
        predat_sp_PSV = flat_hanning_win(tt_sp,predat_sp_PSV,Sp_widewind(1),Sp_widewind(2),3); % 3s taper
        % normalise to unit energy
        normf_sp = predat_sp_PSV(:,1)'*predat_sp_PSV(:,1) + predat_sp_PSV(:,2)'*predat_sp_PSV(:,2);
        predat_sp_PSV = predat_sp_PSV/sqrt(normf_sp);
% -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_S==ir);
        for iir = 1:length(inds)
            final_predata.(Spdat)(inds(iir),1).PSV=predat_sp_PSV; 
            final_predata.(Spdat)(inds(iir),1).tt=tt_sp;
            final_predata.(Spdat)(inds(iir),1).nsamp = length(final_predata.(Spdat)(inds(iir)).PSV);
        end
    end
end % Sp not ccp
end



%% distribute data for different processing (e.g. _lo, _cms)
% for idt = 1:length(par.inv.datatypes)
%     dtype = par.inv.datatypes{idt}; 
%     pdt = parse_dtype( dtype ); 
%     if strcmp(pdt{1},'BW') && (~strcmp(pdt{3},'def') || ~strcmp(pdt{4},'def'))
%         final_predata.(dtype) = final_predata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
%     end
% end


%% ===================  CALCULATE SURPHASE WAVE VELOCITIES  ===================
if any(strcmp(pdtyps(:,1),'SW'))
    % Radial S anis
    xi = zeros(size(final_model.Z));
    xi((final_model.Z>final_model.Zd(1).mu) & (final_model.Z<=final_model.Zd(2).mu)) = final_model.xicrav;
    xi(final_model.Z>final_model.Zd(2).mu) = final_model.ximaav;
    
    modminrun = struct('z',final_model.Z,...
                       'VS',final_model.VSav,...
                        'VP',final_model.VPav,...
                        'rho',final_model.rhoav,...
                        'Sanis',100*(xi-1),...
                        'Panis',zeros(size(final_model.Z)));
                    
    if any(strcmp(pdtyps(:,2),'Ray')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'Ray'),1,'first'));
        par_mineos = struct('R_or_L','R','ID',ID);
        [SW.Ray.phV,SW.Ray.grV] = run_mineos(modminrun,data.(itp{1}).periods,par_mineos,1,0,par.inv.verbose);
    end
    if any(strcmp(pdtyps(:,2),'Lov')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'Lov'),1,'first'));
        par_mineos = struct('R_or_L','L','ID',ID);
        [SW.Lov.phV,SW.Lov.grV] = run_mineos(modminrun,data.(itp{1}).periods,par_mineos,1,0,par.inv.verbose);
    end
    if any(strcmp(pdtyps(:,2),'HV')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'HV'),1,'first'));
        SW.HV.HVr = run_HVkernel(modminrun,data.(itp{1}).periods,'final',1,0,par.inv.verbose);
    end
    for id = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{id}; pdtyp=parse_dtype(dtype); if ~strcmp(pdtyp{1},'SW'),continue; end
        final_predata.(dtype).(pdtyp{3}) = SW.(pdtyp{2}).(pdtyp{3});
    end

end

%% ===================  SP RFs FROM PROPAGATOR MATRIX  ====================
if any(strcmp(pdtyps(strcmp(pdtyps(:,2),'Sp'),3),'ccp'))
    rayp = final_predata.RF_Sp_ccp.rayp;
    % declare samprate of 10 - moot as migrated to depth and interpolated later
    samprate = 10;

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

        % convert to P, SV - use LAYMODEL Vp,Vs to do this precisely
        [predat_sp_PSV(:,1),predat_sp_PSV(:,2)] = ...
            Rotate_XZ_to_PSV(predat_sp_ZRT(:,2),-predat_sp_ZRT(:,1),... % z positive down
            laymodel.Vp(1),laymodel.Vs(1),...
            rayp_sdeg2skm(rayp,laymodel_Suse.zlayb(end)));

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
            tt_sp = []; predat_sp_PSV = nan; % NO GOOD DATA
        end

        %% migrate from time to depth
        %     fprintf(' Migrating RF(z) to RF(t) using %s, rayp = %.3f s/deg\n',par.datprocess.CCP.migratemodel,par.datprocess.CCP.rayp_S)
        simp_final_mod = struct('z',final_model.Z,'VS',final_model.VSav,'VP',final_model.VPav);
        zz_mig = migrate_PorS_conv(tt_sp,simp_final_mod,rayp_sdeg2skm(rayp),0,'Sp');

        % fill in zz_mig "above" surface, corresponding to positive time...
        zz_mig(max(tt_sp)>=tt_sp & tt_sp>0) = -flipud(zz_mig(0>tt_sp & tt_sp>=-max(tt_sp)));
        % discard remaining nans
        kill = isnan(zz_mig);
        zz_mig(kill) = [];
        predat_sp_PSV(kill,:) = [];

        % convert RF_t to RF_z 
        RF_P  = interp1(zz_mig,predat_sp_PSV(:,1),final_predata.RF_Sp_ccp.zz,'linear',0);
        RF_SV  = interp1(zz_mig,predat_sp_PSV(:,2),final_predata.RF_Sp_ccp.zz,'linear',0);

        % taper off daughter = P component
        Zwin = par.datprocess.CCP.Zwin.def;
        taperz = par.datprocess.CCP.taperz;
        RF_Pw = flat_hanning_win(final_predata.RF_Sp_ccp.zz,RF_P,Zwin(1)-taperz/2,Zwin(2)+taperz/2,taperz);    


    else
        final_predata.RF_Sp_ccp.PSV = NaN;    
    end
    % -----------------  PUT INTO DATA STRUCTURE  -----------------
    final_predata.RF_Sp_ccp.PSV = [RF_Pw,RF_SV]; 
    
end

%% ===================  HK-stack: put in values  ====================
if any(strcmp(pdtyps(:,1),'HKstack'))
    HKdat = par.inv.datatypes{find(strcmp(pdtyps(:,1),'HKstack'),1,'first')};
    ik = mindex(data.(HKdat).K,final_model.vpvsav);
    ih = mindex(data.(HKdat).H,zmoh);
    final_predata.(HKdat).H = zmoh;
    final_predata.(HKdat).K = final_model.vpvsav;
    final_predata.(HKdat).E_by_Emax = data.(HKdat).Esum(ik,ih)/maxgrid(data.(HKdat).Esum);
end



end


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


%% ===================  PS RFs FROM PROPAGATOR MATRIX  ====================

if any(strcmp(pdtyps(:,2),'Ps'))
    Psdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Ps'),1,'first')};
    [ unique_rayps_P,irayps_P ] = rayp_vals( [final_predata.(Psdat).rayp] );
    for ir = 1:length(unique_rayps_P)
        rayp = unique_rayps_P(ir);
        samprate = unique([final_predata.(Psdat)(irayps_P==ir).samprate]);
        P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
        [predat_ps,tt_ps] = run_propmat(laymodel,ID,'Ps',samprate, P_inc, par.forc.synthperiod,par.forc.nsamps);
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
    Spdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Sp'),1,'first')};
    [ unique_rayps_S,irayps_S ] = rayp_vals( [final_predata.(Spdat).rayp] );
    for ir = 1:length(unique_rayps_S)
        rayp = unique_rayps_S(ir);
        samprate = unique([final_predata.(Spdat)(irayps_S==ir).samprate]);
        S_inc = rayp2inc(rayp,laymodel.Vs(end),6371-laymodel.zlayb(end));
        [predat_sp,tt_sp] = run_propmat(laymodel,ID,'Sp',samprate, S_inc-0.1, par.forc.synthperiod,par.forc.nsamps);
        predat_sp_ZRT = predat_sp(:,[3,1,2]); % in Z,R,T
        if strcmp(par.forc.PSVorZR,'PSV')
            clear predat_sp_PSV;
            % convert to P, SV
            [predat_sp_PSV(:,1),predat_sp_PSV(:,2)] = ...
                Rotate_XZ_to_PSV(predat_sp_ZRT(:,2),-predat_sp_ZRT(:,1),...
                mean([final_predata.(Spdat).Vp_surf]),mean([final_predata.(Spdat).Vs_surf]),...
                rayp_sdeg2skm(rayp,laymodel.zlayb(end)));
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
    modminrun = struct('z',final_model.Z,...
                       'VS',final_model.VSav,...
                        'VP',final_model.VPav,...
                        'rho',final_model.rhoav,...
                        'Sanis',final_model.Sanis,...
                        'Panis',final_model.Panis);
    
if any(strcmp(pdtyps(:,2),'Ray')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'Ray'),1,'first'));
    [SW.Ray.phV,SW.Ray.grV] = run_mineos(modminrun,data.(itp{1}).periods,'R','final',0,0,par.inv.verbose);
end
if any(strcmp(pdtyps(:,2),'Lov')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'Lov'),1,'first'));
    [SW.Lov.phV,SW.Lov.grV] = run_mineos(modminrun,data.(itp{1}).periods,'L','final',0,0,par.inv.verbose);
end
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id}; pdtyp=parse_dtype(dtype); if ~strcmp(pdtyp{1},'SW'),continue; end
    final_predata.(dtype).phV = SW.(pdtyp{2}).(pdtyp{3});
end

end


end


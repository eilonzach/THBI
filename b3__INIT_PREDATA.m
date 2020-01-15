function [ predata ,laymodel ] = b3__INIT_PREDATA( model,par,trudata,ifplot )
% [ predata ,laymodel ] = b3__INIT_PREDATA( model,par,trudata,ifplot )
% 
%   Wipe trudata to make predata, and layerize model
% 
% INPUTS
%   model   - model structure
%   Kbase   - structure with kernel model, depth kernels, and its phase vels
%   par     - parameters structure
%   predata - data structure with all datatypes 
%   ID      - unique ID for the propmat script to avoid overwriting files
%             if running in parallel.
%   ifplot  - flag with option to plot (1) or not (0)
% 
% OUTPUTS
%   predata - structure identical to input data structure, but with
%             BLANK data, rather than observed data
%%
% An important component is the layerising of the model - conversion of
% continuous model into a bunch of layers, with coarseness partly
% determined by the minimum dVs in any layer (specified as an input). The
% layerised 1D model is also output from this function.

%% ===================  PREPARE DATA STRUCTURE  ===================

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

% [predata.PsRF.Vp_surf] = deal(mean([predata.PsRF.Vp_surf]));
% [predata.PsRF.Vs_surf] = deal(mean([predata.PsRF.Vs_surf]));

%% ===================  ERASE INPUT DATA TO MAKE SURE  ===================
predata = trudata;

indtypes = fieldnames(trudata); % look at all dtypes in predata
for idt = 1:length(indtypes)
    pdt = parse_dtype(indtypes{idt}); % parse data type
    if any(strcmp({'BW','RF'},pdt{1}))
        for irp = 1:length(predata.(indtypes{idt})) % loop over all ray parameters if several
            predata.(indtypes{idt})(irp).PSV(:) = nan; % set to nan
        end
    elseif strcmp({'SW'},pdt{1}) 
        predata.(indtypes{idt}).(pdt{3})(:) = nan; % set to nan
    end
end

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



end


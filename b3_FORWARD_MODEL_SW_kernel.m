function [ predata ] = b3_FORWARD_MODEL_SW_kernel( model,Kbase,par,predata )
%[ predata ] = b3_FORWARD_MODEL_SW_kernel( model,Kbase,par,predata )
% 
%   Do forward model to calculate predicted data.
% 
% INPUTS
%   model   - model structure
%   Kbase   - structure with kernel model, depth kernels, and its phase vels
%   par     - parameters structure
%   data    - obs data structure with all datatypes 
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

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

%% ======================  SURFACE WAVES  ======================
[ modptb ] = calc_Vperturbation(Kbase.modelk,model);

for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    pdtyp = parse_dtype(dtype);
    if ~strcmp(pdtyp{1},'SW'), continue, end
    

%% --------------------  Phase Velocities  --------------------
    switch pdtyp{2}

        case {'Ray','Lov'}

            % calc. perturbation values from 0==>1 and use to calc dc/c

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

%% --------------------  HV ratios  --------------------
   
        case {'HV'}

        % these should be the same each time...
        dvs_vs = linterp(modptb.Z,modptb.dvsav,0.5*(Kbase.(pdtyp{2}).KHV{1}.Z1 + Kbase.(pdtyp{2}).KHV{1}.Z2));
        dvp_vp = linterp(modptb.Z,modptb.dvpav,0.5*(Kbase.(pdtyp{2}).KHV{1}.Z1 + Kbase.(pdtyp{2}).KHV{1}.Z2));
        drh_rh = linterp(modptb.Z,modptb.drho, 0.5*(Kbase.(pdtyp{2}).KHV{1}.Z1 + Kbase.(pdtyp{2}).KHV{1}.Z2));
        zzz = 0.5*(Kbase.(pdtyp{2}).KHV{1}.Z1 + Kbase.(pdtyp{2}).KHV{1}.Z2);

        Np = length(Kbase.(pdtyp{2}).(pdtyp{3}));        
        for ip = 1:Np
            if ~isequal(zzz,0.5*(Kbase.(pdtyp{2}).KHV{ip}.Z1 + Kbase.(pdtyp{2}).KHV{ip}.Z2))
            dvs_vs = linterp(modptb.Z,modptb.dvsav,0.5*(Kbase.(pdtyp{2}).KHV{ip}.Z1 + Kbase.(pdtyp{2}).KHV{ip}.Z2));
            dvp_vp = linterp(modptb.Z,modptb.dvpav,0.5*(Kbase.(pdtyp{2}).KHV{ip}.Z1 + Kbase.(pdtyp{2}).KHV{ip}.Z2));
            drh_rh = linterp(modptb.Z,modptb.drho, 0.5*(Kbase.(pdtyp{2}).KHV{ip}.Z1 + Kbase.(pdtyp{2}).KHV{ip}.Z2));
            end
            dHV = sum(dvs_vs.*Kbase.(pdtyp{2}).KHV{ip}.Kzh_Vs + ...
                      dvp_vp.*Kbase.(pdtyp{2}).KHV{ip}.Kzh_Vp + ...
                      drh_rh.*Kbase.(pdtyp{2}).KHV{ip}.Kzh_rho);
            predata.(dtype).HVr(ip) = Kbase.(pdtyp{2}).(pdtyp{3})(ip) - dHV; % I think the signs are correct here!
        end
 

            
            
    end

end


end


function [Kbase] = make_allkernels(model,Kbase,data,ID,par)
% [Kbase] = make_allkernels(model,Kbase,periods,ID,par)
%   Function to make all surface wave kernels

if isempty(Kbase)
    Kbase = initiate_Kbase;
end

Kbase.modelk = model;
Kbase.itersave=0;
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id}; 
    pdtyp = parse_dtype(dtype); 

    if ~strcmp(pdtyp{1},'SW'), continue; end
    
    if strcmp(pdtyp{2},'HV')
        [HVr_new,HVK_new] = run_HVkernel(model,data.(dtype).periods,ID,1,0,par.inv.verbose);
        Kbase = populate_Kbase( Kbase,dtype,HVr_new,[],{HVK_new} );    
    else
        par_mineos = struct('R_or_L',pdtyp{2},'phV_or_grV',pdtyp{3},'ID',ID);
        [phV,grV,eigfiles] = run_mineos(model,data.(dtype).periods,par_mineos,0,0,par.inv.verbose);
        K = run_kernels(data.(dtype).periods,par_mineos,eigfiles,1,0,par.inv.verbose);
        Kbase = populate_Kbase( Kbase,dtype,phV,grV,{K} );
    end
    
end

end


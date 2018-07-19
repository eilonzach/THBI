function [Kbase] = make_allkernels(model,Kbase,data,ID,par)
% [Kbase] = make_allkernels(model,Kbase,periods,ID,par)
%   Function to make all surface wave kernels

if isempty(Kbase)
    Kbase = initiate_Kbase;
end

Kbase.modelk = model;
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id}; 
    pdtyp = parse_dtype(dtype); 

    if ~strcmp(pdtyp{1},'SW'), continue; end
    
    if strcmp(pdtyp{2},'HV')
        [HVr_new,HVK_new] = run_HVkernel(model,data.(dtype).periods,ID,1,0,par.inv.verbose);
        Kbase = populate_Kbase( Kbase,dtype,HVr_new,[],{HVK_new} );    
    else
        [phV,grV] = run_mineos(model,data.(dtype).periods,pdtyp{2},ID,0,0,par.inv.verbose);
        K = run_kernels(data.(dtype).periods,pdtyp{2},pdtyp{3},ID,1,0,par.inv.verbose);
        Kbase = populate_Kbase( Kbase,dtype,phV,grV,{K} );
    end
    
end

end


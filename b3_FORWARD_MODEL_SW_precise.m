function [ predata,HVK_new,SW ] = b3_FORWARD_MODEL_SW_precise( model,par,predata,ID )
% [ predata,HVK_new,SW ] = b3_FORWARD_MODEL_SW_precise( model,par,predata,ID )
% 
%   Do forward model to calculate predicted data.
% 
% INPUTS
%   model   - model structure
%   Kbase   - structure with kernel model, depth kernels, and its phase vels
%   par     - parameters structure
%   predata - predicted data structure with all datatypes 
%   ID      - unique ID for the propmat script to avoid overwriting files
%             if running in parallel.
% 
% OUTPUTS
%   predata - structure identical to input data structure, but with
%             predicted data, rather than observed data
%%
% An important component is the layerising of the model - conversion of
% continuous model into a bunch of layers, with coarseness partly
% determined by the minimum dVs in any layer (specified as an input). The
% layerised 1D model is also output from this function.

for id = 1:length(par.inv.datatypes)
    allpdytp(id,:)=parse_dtype(par.inv.datatypes{id}); %#ok<AGROW>
end

SW = struct('Ray',struct('phV',[],'grV',[]),'Lov',struct('phV',[],'grV',[]),'HV',struct('HVr',[]));

if any(strcmp(allpdytp(:,2),'Ray')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'Ray'),1,'first'));
    [SW.Ray.phV,SW.Ray.grV] = run_mineos(model,predata.(itp{1}).periods,'R',ID,0,0,par.inv.verbose);
end
if any(strcmp(allpdytp(:,2),'Lov')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'Lov'),1,'first'));
    [SW.Lov.phV,SW.Lov.grV] = run_mineos(model,predata.(itp{1}).periods,'L',ID,0,0,par.inv.verbose);
end
if any(strcmp(allpdytp(:,2),'HV')), itp = par.inv.datatypes(find(strcmp(allpdytp(:,2),'HV'),1,'first'));
    [SW.HV.HVr,HVK_new] = run_HVkernel(model,predata.(itp{1}).periods,['HV_',ID],1,0,par.inv.verbose);
end

for id = 1:length(par.inv.datatypes)

    dtype = par.inv.datatypes{id}; pdtyp=parse_dtype(dtype); 
    if ~strcmp(pdtyp{1},'SW'),continue; end
    switch pdtyp{2}
        case {'Ray','Lov'}
            swk = predata.(dtype).(pdtyp{3}); % record existing grV/phV from kernels
            predata.(dtype).(pdtyp{3}) = SW.(pdtyp{2}).(pdtyp{3});
            swd = predata.(dtype).(pdtyp{3}); % record new, precise grV/phV from mineos
            if par.inv.verbose, unt='m/s'; end
        case 'HV'
            swk = predata.(dtype).HVr; % record existing HVr from kernels
            predata.(dtype).HVr = SW.(pdtyp{2}).(pdtyp{3});
            swd = predata.(dtype).HVr; % record new, precise HVr from Tanimoto script
            if par.inv.verbose, unt=''; end
    end

	if par.inv.verbose
        fprintf('   %s RMS diff is %.4f %s\n',par.inv.datatypes{id},rms(swk-swd),unt); % RMS difference
	end
end


if ifforwardfail(predata,par)
    error('Forward modeling failure thrown during precise SW calculation')
end


end


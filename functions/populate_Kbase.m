function [ Kbase ] = populate_Kbase( Kbase,dtype,phV,grV,K )
%[ Kbase ] = populate_Kbase( Kbase,dtype,phV,grV,K )

pdtyp = parse_dtype(dtype);

if ~iscell(K), K = {K}; end

if ~strcmp(pdtyp{2},'HV') % R or L
    if isempty(Kbase.(pdtyp{2}))
        Kbase.(pdtyp{2}) = struct('phV',phV,'grV',grV,['K',pdtyp{3}(1:2)],K);
    else
        Kbase.(pdtyp{2}).phV = phV;
        Kbase.(pdtyp{2}).grV = grV;
        Kbase.(pdtyp{2}).(['K',pdtyp{3}(1:2)]) = K{1};
    end
elseif strcmp(pdtyp{2},'HV') % HV ratio
    if isempty(Kbase.(pdtyp{2})) 
        Kbase.(pdtyp{2}) = struct('HVr',phV,['K',pdtyp{3}(1:2)],K);
    else
        Kbase.(pdtyp{2}).HVr = phV;
        Kbase.(pdtyp{2}).(['K',pdtyp{3}(1:2)]) = K{1};
    end

end


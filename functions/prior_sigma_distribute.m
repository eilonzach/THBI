function par = prior_sigma_distribute(par,trudata)
%par = prior_sigma_distribute(par,trudata)
%   set prior sigma as geometric mean of data sigma

for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt}; 
    pdt = parse_dtype( dtype ); 
    if isfield(trudata.(dtype),'sigma') && ~isnan(geomean(trudata.(dtype).sigma))
        par.mod.data.prior_sigma.(pdt{1}).(pdt{2}).(pdt{3}) = geomean(trudata.(dtype).sigma);
    end
end

end


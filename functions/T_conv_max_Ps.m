function tmax = T_conv_max_Ps(par)
%Tmax = T_conv_max_Ps(par)
%   Function to find the max P-s conversion time included in the time
%   window used for this inversion;

tmax = 0;
for id = 1:length(par.inv.datatypes)
    dtype = parse_dtype(par.inv.datatypes{id}); 
    if strcmp(dtype{1},'BW') && strcmp(dtype{2},'Ps')
        tmax = max([tmax,par.datprocess.(dtype{2}).Twin.(dtype{3})(2)]);
    end
end


end


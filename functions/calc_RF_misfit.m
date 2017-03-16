function [ misfit ] = calc_misfit( dobs,gm,invCe )
%[ misfit ] = calc_misfit( dobs,gm,invCe )
%   Detailed explanation goes here

misfit = (dobs-gm)'*invCe*(dobs-gm);

end


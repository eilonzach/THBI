function [sub_gradient,supra_gradient] = moho_adj_gradients(model,dz)
% [sub_gradient,supra_gradient] = moho_adj_gradients(model)
%   
%  calculate the gradient of the velocity profile immediately above the
%  moho and immediately beneath it. Positive values indicate increase in
%  velocity with depth

% depth over which to compute gradient
if nargin < 2 || isempty(dz)
    dz = 2*diff(model.z(1:2)); 
end

% sub-moho gradient
z_sub = find((model.z >= model.zmoh) & (model.z <= model.zmoh+dz));
z_sub(1) = []; % get rid of top of moho
sub_gradient = mean(diff(model.VS(z_sub))./diff(model.z(z_sub)));

% supra-moho gradient
z_sup = find((model.z >= model.zmoh - dz) & (model.z <= model.zmoh));
z_sup(end) = []; % get rid of bottom of moho
supra_gradient = mean(diff(model.VS(z_sup))./diff(model.z(z_sup)));


end


function [ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels_collated,par )
% [ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels,par,ifsave )

if nargin < 3 || isempty(ifsave)
    ifsave = 0;
end

am = allmodels_collated;
if ~isfield(am,'chain'), am = dealto(am,'chain',1); end

% don't use pre burn-in models! Just in case
am([am.bestmods]'==false) =  [];

fprintf('  > Building suite of models\n')

%% =====================  Suite of models  =====================
Zsuite = sort(unique([[0:2:par.mod.maxz]';[-15:0.2:15]'+mean([am.zsed]);[-15:0.2:15]'+mean([am.zmoh])]));
Zsuite(Zsuite<0) = [];
Nz = length(Zsuite);
Nm = length(am);
fprintf('  > Resolving suite of models onto common basis \n')
VSsuite = zeros(Nz,Nm);
VPsuite = zeros(Nz,Nm);
rhosuite = zeros(Nz,Nm);

for ii = 1:Nm
VSsuite(:,ii) = linterp(am(ii).z,am(ii).VS,Zsuite);
VPsuite(:,ii) = linterp(am(ii).z,am(ii).VP,Zsuite);
rhosuite(:,ii) = linterp(am(ii).z,am(ii).rho,Zsuite);
end

suite_of_models.Z = Zsuite;
suite_of_models.VS = VSsuite;
suite_of_models.VP = VPsuite;
suite_of_models.rho = rhosuite;
suite_of_models.chain = [am.chain];




end









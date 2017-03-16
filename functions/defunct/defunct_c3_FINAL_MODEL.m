function [ final_model,suite_of_models ] = c3_FINAL_MODEL( posterior,allmodels,par,ifsave )
%C3_FINAL_MODEL Summary of this function goes here
%   Detailed explanation goes here


%% gather models at each depth
Zd1 = [-5:0.2:5]+mean(posterior.zsed);
Zd2 = [-5:0.2:5]+mean(posterior.zmoh);
Zgather = unique(sort([[0:1:par.mod.maxz],Zd1,Zd2])'); Zgather(Zgather<0)=[];
Nz = length(Zgather);

fprintf('  > Resolving all models onto common basis \n')
VSgather = zeros(Nz,posterior.Nstored);
VPgather = zeros(Nz,posterior.Nstored);
rhogather = zeros(Nz,posterior.Nstored);
for ii = 1:posterior.Nstored
VSgather(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VS,Zgather);
VPgather(:,ii) = linterp(allmodels(ii).z,allmodels(ii).VP,Zgather);
rhogather(:,ii) = linterp(allmodels(ii).z,allmodels(ii).rho,Zgather);
end

fprintf('  > Fitting parms at each depth \n')
VSbest = zeros(Nz,1);
VPbest = zeros(Nz,1);
rhobest = zeros(Nz,1);
VSstd = zeros(Nz,1);
VPstd = zeros(Nz,1);
rhostd = zeros(Nz,1);
for iz = 1:Nz
    X = midpts([0:0.1:10]');
    Nvs = hist(VSgather(iz,:),X);
    [VSstd(iz), VSbest(iz)] = gaussfit( X, Nvs );
    Nvp = hist(VPgather(iz,:),X);
    [VPstd(iz), VPbest(iz)] = gaussfit( X, Nvp );
    Nrho = hist(rhogather(iz,:),X);
    [rhostd(iz), rhobest(iz)] = gaussfit( X, Nrho );
end

final_model.Zgather = Zgather;
final_model.VSbest = VSbest;
final_model.VSstd = VSstd;
final_model.VPbest = VPbest;
final_model.VPstd = VPstd;
final_model.rhobest = rhobest;
final_model.rhostd = rhostd;

suite_of_models.Z = Zgather;
suite_of_models.VS = VSgather;
suite_of_models.VP = VPgather;
suite_of_models.rho = rhogather;

if ifsave
    save('final_model','final_model')
end




end









function [ final_model,suite_of_models ] = c4_FINAL_MODEL( posterior,allmodels,par,ifsave,ofile )
% [ final_model ] = c4_FINAL_MODEL( posterior,allmodels,par,ifsave )
%   Detailed explanation goes here

if nargin < 4 || isempty(ifsave)
    ifsave = 0;
end

if nargin < 5 || isempty(ofile)
    ofile = './final_model';
end
 
o4 = ones(1,4);
 
% don't use pre burn-in models!
allmodels([allmodels.bestmods]'==false) =  []; 


%% ===================== GATHER MODELS, FIT WITH GAUSSIANS ================
%% DISCONTINUITIES
Zd = struct('mu',[],'std',[]); 
% crust/sed
X = midpts([0:0.1:par.mod.sed.hmax]);  XX = [0*o4,X,par.mod.sed.hmax*o4];
Nds = hist(posterior.zsed,X); NNds = [Nds(1)*o4,Nds,Nds(end)*o4];
if sum(Nds~=0)==1
    Zd(1).mu = mean(posterior.zsed); Zd(1).std = par.mod.dz/2;
else
    [Zd(1).std,Zd(1).mu] = gaussfit( XX, NNds, (X(end)-X(1))/4,X(Nds==max(Nds))  );
end
if Zd(1).mu < 0, Zd(1).mu = 0; end

% Moho
X = midpts([0:0.5:par.mod.sed.hmax+par.mod.crust.hmax]);
Nds = hist(posterior.zmoh,X);
if sum(Nds~=0)==1
    Zd(2).mu = mean(posterior.zmoh); Zd(2).std = par.mod.dz/2;
else
    [Zd(2).std,Zd(2).mu] = gaussfit( X, Nds,(X(end)-X(1))/4, mean(X(Nds==max(Nds))) );
end

%% Loop through depths obtaining max/min
Z = sort(unique([[0:par.mod.dz:par.mod.maxz],...
                 Zd(1).mu+[0,-3*Zd(1).std:0.1:3*Zd(1).std],...
                 Zd(2).mu+[0,-3*Zd(2).std:0.2:3*Zd(2).std]]')); 
             Z(Z<0)=[];
Nz = length(Z);
Nm = posterior.Nstored;

fprintf('  > Gathering models (layer-wise)\n    & fitting gaussians to each parm/depth\n      ')

variable = {'VS','VP','rho'};
for iv = 1:length(variable)
    fprintf('%s... ',variable{iv});

varz = zeros(Nz,Nm);
for ii = 1:Nm
    varz(:,ii) = linterp(allmodels(ii).z,allmodels(ii).(variable{iv}),Z);
end 

lb = zeros(Nz,1); ub = lb; 
l1std = lb;  u1std = lb; 
l2std = lb;  u2std = lb; 
mdn = lb;

for iz = 1:Nz
    vord = sort(varz(iz,:));
    lb(iz) = vord(round(0.005*Nm));     % lower bound (not quite min, to avoid wacky ones)
    ub(iz) = vord(round(0.995*Nm));     % upper bound (not quite max, to avoid wacky ones)
    l1std(iz) = vord(round(0.3173*Nm)); % 1 std below median
    u1std(iz) = vord(round(0.6827*Nm)); % 1 std above median
    l2std(iz) = vord(round(0.0455*Nm)); % 2 std below median
    u2std(iz) = vord(round(0.9545*Nm)); % 2 std above median
    mdn(iz) = vord(round(0.5*Nm));      % median value
end

final_model.Z = Z;
final_model.([variable{iv},'av']) = mdn;
final_model.([variable{iv},'minmax']) = [lb ub];
final_model.([variable{iv},'sig1']) = [l1std u1std];
final_model.([variable{iv},'sig2']) = [l2std u2std];

suite_of_models.(variable{iv}) = varz;

end % loop on variables
fprintf('\n');


%%  HYPERPARAMETERS 
fns = fieldnames(posterior.datahparm);
hyperparms = struct([]);
for ihp = 1:length(fns)
    dtype = fns{ihp};
    for it = 1:size(posterior.datahparm.(dtype),2)
        postvar = lognfit(posterior.datahparm.(dtype)(:,it),0.05);
        hyperparms(1).(dtype).mu_log10(it,1) = postvar(1)/log(10);
        hyperparms(1).(dtype).std_log10(it,1) = postvar(2)/log(10);
    end
end

%% ===================== add other parameters ================

final_model.hyperparms = hyperparms;
final_model.Zd = Zd;


if ifsave
    save(ofile,'final_model')
end




end









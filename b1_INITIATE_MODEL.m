function model = b1_INITIATE_MODEL(par,selev,ifplot)
% model = b1_INITIATE_MODEL(par)
% 
% Function to make the nstas x 1-D model as per the parameters defined in
% the input structure par. If relevant, give station elevation in km.
% 
% % model parameters:
%   1. hsed
%   2. vs_sed_top
%   3. vs_sed_bottom
%   4. hcrust
%   5. vpvs_crust (may be fixed)
%   ++ V_knots_crust (there are N_knots_crust of these)
%   ++ V_knots_mantle (there are N_knots_mantle of these)
%
%  ==> 5 + Nkc + Nkm parms 
% 
%       (c.f. 13 for Shen and Ritz, who fixed vpvs and Nkc=4, Nkm=5

if nargin <2 || isempty(selev)
    selev = zeros(par.mod.nstas,1);
end
if nargin <3 || isempty(ifplot)
    ifplot = 0;
end

mod = par.mod;

%% resolve important values from prior distributions
h_sed = random('unif',mod.sed.hmin,mod.sed.hmax);
h_crust = random('unif',mod.crust.hmin,mod.crust.hmax);

k_crust = mod.crust.kmin + ceil(random('unid',mod.crust.kmax-mod.crust.kmin+1)/2)-1;
k_mantle = mod.mantle.kmin + ceil(random('unid',mod.mantle.kmax-mod.mantle.kmin+1)/3)-1;

vs_sed = random('unif',mod.sed.vsmin,mod.sed.vsmax,2,1);
kvs_crust = random('unif',mod.crust.vsmin,mod.crust.vsmax,k_crust+1,1);
kvs_mantle = random('unif',mod.mantle.vsmin,mod.mantle.vsmax,k_mantle+1,1);

vpvs_crust = random('unif',mod.crust.vpvsmin,mod.crust.vpvsmax);

xi_crust = random('unif',mod.crust.ximin,mod.crust.ximax);
xi_mantle = random('unif',mod.mantle.ximin,mod.mantle.ximax);


% impose monotonic increase conditions:
vs_sed = sort(vs_sed);
% kvs_crust = sort(kvs_crust);

%% SEDIMENT PARMS
sed = struct('h',h_sed,'VS',vs_sed);

%% CRUST PARMS
cminz = h_sed;
cmaxz = h_sed+h_crust;
zc = unique([cminz:mod.dz:cmaxz,cmaxz])';

% set up splines
fcknots = linspace(0,1,k_crust)';
cknots = cminz + (cmaxz-cminz)*fcknots; % linearly spaced knots

[spbasis]=make_splines(cknots,par,zc,zc);

crust = struct('h',h_crust,'Nsp',k_crust+1,'Nkn',k_crust,'VS_sp',kvs_crust,'vpvs',vpvs_crust,'xi',xi_crust,'splines',spbasis,'knots',cknots,'fknots',fcknots,'z_sp',zc);


%% MANTLE PARMS
mminz = h_sed+h_crust;
mmaxz = mod.maxz + selev;
zm = unique([mminz:mod.dz:mmaxz,mmaxz])';

% set up splines
fmknots_ = linspace(0,1,k_mantle-1)';  % linearly spaced knots
mknots = [mminz + (mod.maxkz-mminz)*fmknots_ ; mmaxz];
fmknots = (mknots - mminz)./(mmaxz-mminz);

[spbasis] = make_splines(mknots,par,zm,zm);

mantle = struct('Nsp',k_mantle+1,'Nkn',k_mantle,'VS_sp',kvs_mantle,'xi',xi_mantle,'splines',spbasis,'knots',mknots,'fknots',fmknots,'z_sp',zm);

%% DATA PARMS
hparm = struct([]);
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    [ ~,sigma_prior ] = get_sigma_min_prior( dtype,par );
    hparm(1).(['sig_',dtype]) = sigma_prior;
end

%% OVERALL PARMS
M = 2 + 1 + k_crust + k_mantle ...
      + double(mod.sed.hmin~=mod.sed.hmax) ...
      + double(mod.crust.vpvsmin~=mod.crust.vpvsmax) ...
      + double(mod.crust.ximin~=mod.crust.ximax) ...
      + double(mod.mantle.ximin~=mod.mantle.ximax);

%% MODEL WITH ALL PARMS
model = struct('sedmparm',sed,'crustmparm',crust,'mantmparm',mantle,'datahparm',hparm,...
               'M',M,'selev',selev);
           
%% TURN PARMS INTO REAL MODEL
model = make_mod_from_parms(model,par);

% fknots = linspace(0,1,model.mantmparm.Nsp-1);
% [ aa,bb,cc ] = make_splines_fraction( model,par,fknots,'mantle',model.z(15:end),model.VS(15:end) )
% [ dd,ee ] = make_splines( cc,par,model.z(15:end),model.VS(15:end) )


%% plot
if ifplot
    figure(31),clf
    plot(model.VS,model.z,'b',model.VP,model.z,'r',model.rho,model.z,'k')
    set(gca,'ydir','reverse')
end
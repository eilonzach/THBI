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
%   5. N_knots_crust (may be fixed)
%   ++. V_knots_crust (there are N_knots_crust of these)
%   6. N_knots_mantle (may be fixed)
%   ++. V_knots_mantle (there are N_knots_mantle of these)
%   7. vpvs_mantle
%
%  ==> 7 + Nkc + Nkm parms (fewer if N_knots are fixed)
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

k_crust = mod.crust.kmin + random('unid',mod.crust.kmax-mod.crust.kmin+1)-1;
k_mantle = mod.mantle.kmin + random('unid',mod.mantle.kmax-mod.mantle.kmin+1)-1;

vs_sed = random('unif',mod.sed.vsmin,mod.sed.vsmax,2,1);
kvs_crust = random('unif',mod.crust.vsmin,mod.crust.vsmax,k_crust,1);
kvs_mantle = random('unif',mod.mantle.vsmin,mod.mantle.vsmax,k_mantle,1);

vpvs_mantle = random('unif',mod.mantle.vpvsmin,mod.mantle.vpvsmax);
vpvs_crust = vpvs_mantle;

% impose monotonic increase conditions:
vs_sed = sort(vs_sed);
kvs_crust = sort(kvs_crust);

%% SEDIMENTS
zs = [0,h_sed]';
vp_sed = sed_vs2vp(vs_sed);
rho_sed = sed_vs2rho(vs_sed);
sed = struct('h',h_sed,'vstop',vs_sed(1),'vsbot',vs_sed(2));

%% CRUST
cminz = h_sed;
cmaxz = h_sed+h_crust;
zc = unique([cminz:mod.dz:cmaxz,cmaxz])';

% set up splines
dzsp = (cmaxz-cminz)/(k_crust-2);
knots = [repmat(cminz,1,3),cminz:dzsp:cmaxz,repmat(cmaxz,1,3)]';
sp = fastBSpline.lsqspline(knots,2,linspace(cminz,cmaxz,k_crust)',kvs_crust); % dummy velocities as placeholder
spbasis = sp.getBasis(zc); spbasis = spbasis(:,2:end-1);

vs_crust = sum(spbasis*diag(kvs_crust),2);
vp_crust = vs_crust*vpvs_crust;
rho_crust = sed_vs2rho(vs_crust); % use same expression as for sed (assume not ultramafic or other poorly-fit rocks)

crust = struct('h',h_crust,'Nkn',k_crust,'VS_kn',kvs_crust,'splines',spbasis);

%% MANTLE
mminz = h_sed+h_crust;
mmaxz = mod.maxz + selev;
zm = unique([mminz:mod.dz:mmaxz,mmaxz])';

% set up splines
dzsp = (mmaxz-mminz)/(k_mantle-2);
knots = [repmat(mminz,1,3),mminz:dzsp:mmaxz,repmat(mmaxz,1,3)]';
sp = fastBSpline.lsqspline(knots,2,linspace(mminz,mmaxz,k_mantle)',kvs_mantle); % dummy velocities as placeholder
spbasis = sp.getBasis(zm); spbasis = spbasis(:,2:end-1);

vs_mantle = sum(spbasis*diag(kvs_mantle),2);
vp_mantle = vs_mantle*vpvs_mantle;
rho_mantle = mantle_vs2rho(vs_mantle,zm );

mantle = struct('Nkn',k_mantle,'VS_kn',kvs_crust,'splines',spbasis,'vpvs',vpvs_mantle);


%% COLLATE
zz = [zs;zc;zm];
zz0 = zz-selev; % true depth, from sea level/ref. ellipsoid
Nz = length(zz);
zsed = h_sed;
zmoh = h_sed+h_crust;
vs = [vs_sed;vs_crust;vs_mantle];
vp = [vp_sed;vp_crust;vp_mantle];
rho = [rho_sed;rho_crust;rho_mantle];

%% OUTPUT
model = struct('z',zz,'z0',zz0,'VS',vs,'VP',vp,'rho',rho,...
               'Nz',Nz,'zsed',zsed,'zmoh',zmoh,'selev',selev,...
               'fdVSsed',100*diff(vs(zz==zsed))./mean(vs(zz==zsed)),...
               'fdVSmoh',100*diff(vs(zz==zmoh))./mean(vs(zz==zmoh)),...
               'Sanis',zeros(Nz,1),'Panis',zeros(Nz,1),...
               'sedmparm',sed,'crustmparm',crust,'mantmparm',mantle,...
               'M',7+k_crust+k_mantle);



%% plot
if ifplot
    figure(31),clf
    plot(vs,zz,'b',vp,zz,'r',rho,zz,'k')
    set(gca,'ydir','reverse')
end
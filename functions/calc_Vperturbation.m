function [ modptb ] = calc_Vperturbation( model0,model1,ifplot )
%[ modptb ] = calc_Vperturbation( model0,model1,ifplot )
% 
% Function to calculate the fractional perturbation of model1 from the
% reference model0, as a function of depth, for vsv,vsh,vpv,vph,rho
% 
% perturbations are returned for all depths in corresponding card file
% (from surface to Earth centre)

if nargin <3 || isempty(ifplot)
    ifplot=false;
end

%% Radial anisotropy
xi0  = 1 + model0.Sanis/100;  % assumes Sanis is a percentage of anis about zero
phi0 = 1 + model0.Panis/100;  % assumes Panis is a percentage of anis about zero
[ vsv0,vsh0 ] = VsvVsh_from_VsXi( model0.VS,xi0 );
[ vpv0,vph0 ] = VpvVph_from_VpPhi( model0.VP,phi0 );

xi1  = 1 + model1.Sanis/100;  % assumes Sanis is a percentage of anis about zero
phi1 = 1 + model1.Panis/100;  % assumes Panis is a percentage of anis about zero
[ vsv1,vsh1 ] = VsvVsh_from_VsXi( model1.VS,xi1 );
[ vpv1,vph1 ] = VpvVph_from_VpPhi( model1.VP,phi1 );


%% get in card format inc. all depths in PREM
card0 = write_cardfile([],model0.z,vpv0,vsv0,model0.rho,[],[],vph0,vsh0);
card1 = write_cardfile([],model1.z,vpv1,vsv1,model1.rho,[],[],vph1,vsh1);

%% parse depths and nodes
zz0 = card0.depth;
zz1 = card1.depth;
N = length(zz0); % this many kernel rows

% find coincident nodes
[zzbotha,zbothi0a,zbothi1a] = intersect(zz0,zz1,'stable'); % gets first indices of discs
[zzbothb,zbothi0b,zbothi1b] = intersect(zz0,zz1,'legacy'); % gets second indices of discs
zbothi0 = union(zbothi0a,zbothi0b);
zbothi1 = union(zbothi1a,zbothi1b);

zdisc0 = card0.depth(diff(card0.depth)==0);
zdisc1 = card1.depth(diff(card1.depth)==0);

% account for 1's discontinuities matching single nodes in 0.
if length(zbothi0)~=length(zbothi1)
    d1 = setdiff(zdisc1,zdisc0);
    if any([d1])
        for ii = 1:length(d1)
            zbothi0 = sort([zbothi0;find(zz0==d1(ii))]);
        end
    end
end
if ~isequal(zz0(zbothi0),zz1(zbothi1))
    error('something wrong with indices'); 
end

%find non-coincident points
zdiffi0 = setdiff([1:N]',zbothi0);


%% calc dvals
% initialise dvals
dvsv = nan(N,1);
dvsh = nan(N,1);
dvpv = nan(N,1);
dvph = nan(N,1);
drho = nan(N,1);

% insert dvals for coincident points
dvsv(zbothi0) = card1.vsv(zbothi1)./card0.vsv(zbothi0) - 1;
dvsh(zbothi0) = card1.vsh(zbothi1)./card0.vsh(zbothi0) - 1;
dvpv(zbothi0) = card1.vpv(zbothi1)./card0.vpv(zbothi0) - 1;
dvph(zbothi0) = card1.vph(zbothi1)./card0.vph(zbothi0) - 1;
drho(zbothi0) = card1.rho(zbothi1)./card0.rho(zbothi0) - 1;

% insert dvals for non-coincident points (resolve onto card0 basis)
dvsv(zdiffi0) = linterp(zz1,card1.vsv,zz0(zdiffi0))./card0.vsv(zdiffi0) - 1;
dvsh(zdiffi0) = linterp(zz1,card1.vsh,zz0(zdiffi0))./card0.vsh(zdiffi0) - 1;
dvpv(zdiffi0) = linterp(zz1,card1.vpv,zz0(zdiffi0))./card0.vpv(zdiffi0) - 1;
dvph(zdiffi0) = linterp(zz1,card1.vph,zz0(zdiffi0))./card0.vph(zdiffi0) - 1;
drho(zdiffi0) = linterp(zz1,card1.rho,zz0(zdiffi0))./card0.rho(zdiffi0) - 1;

% fix nans
dvsv(card0.vsv==0) = 0; if any(isnan(dvsv)), error('nans in dvsv'); end
dvsh(card0.vsh==0) = 0; if any(isnan(dvsh)), error('nans in dvsh'); end
dvpv(card0.vpv==0) = 0; if any(isnan(dvpv)), error('nans in dvpv'); end
dvph(card0.vph==0) = 0; if any(isnan(dvph)), error('nans in dvph'); end
drho(card0.rho==0) = 0; if any(isnan(drho)), error('nans in drho'); end


%% output
modptb =struct('Z',zz0,'dvsv',dvsv,'dvsh',dvsh,'dvpv',dvpv,'dvph',dvph,'drho',drho);

%% plot
if ifplot
    figure(44); clf; set(gcf,'pos',[331 384 1537 614])
    ax1 = subplot(1,5,1); 
    ax2 = subplot(1,5,2); 
    ax3 = subplot(1,5,3); 
    ax4 = subplot(1,5,4); 
    ax5 = subplot(1,5,5); 
    
    plot(ax1,dvsv,zz0,'linewidth',2);
    plot(ax2,dvsh,zz0,'linewidth',2);
    plot(ax3,dvpv,zz0,'linewidth',2);
    plot(ax4,dvph,zz0,'linewidth',2);
    plot(ax5,drho,zz0,'linewidth',2);
    
    set(ax1,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax2,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax3,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax4,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax5,'ydir','reverse','fontsize',15,'ylim',[0 500])
    
    title(ax1,'Vsv','fontsize',20);xlabel(ax1,'$\mathbf{dV_{SV}/V_{SV}}$','interpreter','latex','fontsize',19)
    title(ax2,'Vsh','fontsize',20);xlabel(ax2,'$\mathbf{dV_{SH}/V_{SH}}$','interpreter','latex','fontsize',19)
    title(ax3,'Vpv','fontsize',20);xlabel(ax3,'$\mathbf{dV_{PV}/V_{PV}}$','interpreter','latex','fontsize',19)
    title(ax4,'Vph','fontsize',20);xlabel(ax4,'$\mathbf{dV_{PH}/V_{PH}}$','interpreter','latex','fontsize',19)
    title(ax5,'rho','fontsize',20);xlabel(ax5,'$\mathbf{d\rho/\rho}$','interpreter','latex','fontsize',19)
    ylabel(ax1,'\textbf{Depth (km)}','interpreter','latex','fontsize',19)
end

end


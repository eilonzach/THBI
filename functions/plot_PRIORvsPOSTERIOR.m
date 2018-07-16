function plot_PRIORvsPOSTERIOR(prior,posterior,par,ifsave,ofile)
% plot_PRIORvsPOSTERIOR(prior,posterior,par,ifsave,ofile)
%   
% function to compare prior and posterior models w/ some relevant metrics

if nargin<4 || isempty(ifsave)
    ifsave = 0; % default is not to save
end
if nargin<5 || isempty(ofile)
    ofile = 'figs/priorvsposterior_fig.pdf';
end
 
    
%% plot prior
figure(91), clf, set(gcf,'pos',[43 141 1774 961]);
cls = get(groot,'defaultAxesColorOrder');

% try % this will fail if no sed variation because it is fixed
% subplot(341), cla, hold on
% X = midpts(linspace(par.mod.sed.hmin,par.mod.sed.hmax,20));
% No = hist(posterior.zsed,X)/posterior.Nstored;
% Ni = hist(prior.zsed,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14), title('Sed thickness (km)','fontsize',16)
% 
% subplot(342), cla, hold on
% X = midpts(linspace(par.mod.sed.vsmin,par.mod.sed.vsmax,20));
% No = hist(posterior.VSsedtop,X)/posterior.Nstored;
% Ni = hist(prior.VSsedtop,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14), title('Vs seds top (km/s)','fontsize',16)
% 
% subplot(343), cla, hold on
% X = midpts(linspace(par.mod.sed.vsmin,par.mod.sed.vsmax,20));
% No = hist(posterior.VSsedbot,X)/posterior.Nstored;
% Ni = hist(prior.VSsedbot,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14), title('Vs seds bot (km/s)','fontsize',16)
% 
% subplot(344), cla, hold on
% title('fractional dVs at sed/crust (%)','fontsize',16)
% X = midpts(linspace(0,150,20));
% No = hist(posterior.fdVSsed,X)/posterior.Nstored;
% Ni = hist(prior.fdVSsed,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14)
% catch
%     delete(subplot(341));delete(subplot(342));delete(subplot(343));delete(subplot(344));
% end
% 
% subplot(345), cla, hold on
% title('Moho depth (km)','fontsize',16)
% X = midpts(linspace(par.mod.sed.hmin+par.mod.crust.hmin,par.mod.sed.hmax+par.mod.crust.hmax,20));
% No = hist(posterior.zmoh,X)/posterior.Nstored;
% Ni = hist(prior.zmoh,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14,'xlim',[par.mod.sed.hmin+par.mod.crust.hmin par.mod.sed.hmax+par.mod.crust.hmax])
% 
% subplot(346), cla, hold on
% X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
% No = hist(posterior.VScrusttop,X)/posterior.Nstored;
% Ni = hist(prior.VScrusttop,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14,'xlim',[par.mod.crust.vsmin,par.mod.crust.vsmax]), title('Vs crust top (km/s)','fontsize',16)
% 
% subplot(347), cla, hold on
% title('Vs crust bot (km/s)','fontsize',16)
% X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
% No = hist(posterior.VScrustbot,X)/posterior.Nstored;
% Ni = hist(prior.VScrustbot,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14)
% 
% subplot(348), cla, hold on
% X = midpts(linspace(0,30,20));
% No = hist(posterior.fdVSmoh,X)/posterior.Nstored;
% Ni = hist(prior.fdVSmoh,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14), title('fractional dVs at moho (%)','fontsize',16)
% 
% subplot(349), cla, hold on
% title('crust Vp/Vs ratio','fontsize',16)
% X = midpts(linspace(par.mod.crust.vpvsmin,par.mod.crust.vpvsmax,20));
% No = hist(posterior.vpvs,X)/posterior.Nstored;
% Ni = hist(prior.vpvs,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14,'xlim',[par.mod.crust.vpvsmin,par.mod.crust.vpvsmax]) 
% 
% 
% %% velocity at depth
% X = midpts(linspace(par.mod.mantle.vsmin,par.mod.mantle.vsmax,20));
% zdo = [50,65,75,90,120,150]; for iz = 1:length(zdo), izdo(iz) = crossing(posterior.zmantle,[],zdo(iz));end
% for iz = 1:6
% subplot(3,8,18+iz), cla, hold on
% No = hist(posterior.VSmantle(:,izdo(iz)),X)/posterior.Nstored;
% Ni = hist(prior.VSmantle(:,izdo(iz)),X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% % legend(num2str(model_summary.zmantle(:)),'location','northwest')
% set(gca,'fontsize',14,'xlim',[3.5 4.9]), title(sprintf('Vs at %.0f km',prior.zmantle(izdo(iz))),'fontsize',16)
% end

%% Set up plot
figure(19), clf, set(gcf,'pos',[43 141 1774 661]);
cls = get(groot,'defaultAxesColorOrder');

% subplots in grid
% Nx = 4; dx = (0.9-0.03*(Nx-1))/Nx;
% Ny = 2; dy = (0.84-0.07*(Ny-1))/Ny;
% for iy = 1:2, for ix = 1:4;
% ax(iy,ix) = axes(gcf,'pos',[(0.07 + (dx+0.03)*(ix-1)) (0.07 + (dy+0.07)*(iy-1)) dx dy]); hold on, 
% end,end

% plot each data type
subplot(241), cla, hold on
title('Moho depth (km)','fontsize',16)
X = midpts(linspace(par.mod.sed.hmin+par.mod.crust.hmin,par.mod.sed.hmax+par.mod.crust.hmax,20));
No = hist(posterior.zmoh,X)/posterior.Nstored;
Ni = hist(prior.zmoh,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.sed.hmin+par.mod.crust.hmin par.mod.sed.hmax+par.mod.crust.hmax])

subplot(242), cla, hold on
X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
No = hist(posterior.VScrusttop,X)/posterior.Nstored;
Ni = hist(prior.VScrusttop,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.crust.vsmin,par.mod.crust.vsmax]), title('Vs crust top (km/s)','fontsize',16)

subplot(243), cla, hold on
title('Vs crust bot (km/s)','fontsize',16)
X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
No = hist(posterior.VScrustbot,X)/posterior.Nstored;
Ni = hist(prior.VScrustbot,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14)

subplot(244), cla, hold on
X = midpts(linspace(0,30,20));
No = hist(posterior.fdVSmoh,X)/posterior.Nstored;
Ni = hist(prior.fdVSmoh,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14), title('fractional dVs at moho (%)','fontsize',16)

subplot(245), cla, hold on
title('crust Vp/Vs ratio','fontsize',16)
X = midpts(linspace(par.mod.crust.vpvsmin,par.mod.crust.vpvsmax,20));
No = hist(posterior.vpvs,X)/posterior.Nstored;
Ni = hist(prior.vpvs,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'xlim',[par.mod.crust.vpvsmin,par.mod.crust.vpvsmax]) 


%% velocity at depth
X = midpts(linspace(par.mod.crust.vsmin,par.mod.mantle.vsmax,30));
zdo = [53,65,75,90,120,150]; for iz = 1:length(zdo), izdo(iz) = crossing(posterior.zatdep,[],zdo(iz));end
for iz = 1:6
subplot(2,8,10+iz), cla, hold on
No = hist(posterior.VSmantle(:,izdo(iz)),X)/posterior.Nstored;
Ni = hist(prior.VSmantle(:,izdo(iz)),X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% legend(num2str(model_summary.zmantle(:)),'location','northwest')
set(gca,'fontsize',14,'xlim',[3.5 4.9]), title(sprintf('Vs at %.0f km',prior.zatdep(izdo(iz))),'fontsize',16)
end

%% title
htit = title_custom([par.sta,' ',par.nwk],0.05,'fontweight','bold','fontsize',25);

if ifsave
    save2pdf(19,ofile,'/');
end


end


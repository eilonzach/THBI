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
 
%% Set up plot
figure(19), clf, set(gcf,'pos',[43 141 1774 961]);
cls = get(groot,'defaultAxesColorOrder');

%% Moho depth
subplot(341), cla, hold on
X = midpts(linspace(par.mod.sed.hmin+par.mod.crust.hmin,par.mod.sed.hmax+par.mod.crust.hmax,20));
No = hist(posterior.zmoh,X)/posterior.Nstored;
Ni = hist(prior.zmoh,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'xlim',[par.mod.sed.hmin+par.mod.crust.hmin par.mod.sed.hmax+par.mod.crust.hmax],'ylim',[0 axlim(gca,4)])
title('Moho depth (km)','fontsize',16)

%% V crust bottom (commented out top)
% subplot(342), cla, hold on
% X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
% No = hist(posterior.VScrusttop,X)/posterior.Nstored;
% Ni = hist(prior.VScrusttop,X)/prior.Nstored;
% bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
% bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
%         'xlim',[par.mod.crust.vsmin,par.mod.crust.vsmax],'ylim',[0 axlim(gca,4)])
% title('Vs crust top (km/s)','fontsize',16)

subplot(342), cla, hold on
X = midpts(linspace(par.mod.crust.vsmin,par.mod.crust.vsmax,20));
No = hist(posterior.VScrustbot,X)/posterior.Nstored;
Ni = hist(prior.VScrustbot,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'xlim',[par.mod.crust.vsmin,par.mod.crust.vsmax],'ylim',[0 axlim(gca,4)])
title('Vs crust bot (km/s)','fontsize',16)

%% Moho dVs
subplot(343), cla, hold on
X = midpts(linspace(0,30,20));
No = hist(posterior.fdVSmoh,X)/posterior.Nstored;
Ni = hist(prior.fdVSmoh,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'ylim',[0 axlim(gca,4)])
title('fractional dVs at Moho (%)','fontsize',16)

%% H-K value
subplot(344), cla, hold on
% Xh = linspace(par.mod.sed.hmin+par.mod.crust.hmin,par.mod.sed.hmax+par.mod.crust.hmax,30);
% Xk = linspace(par.mod.crust.vpvsmin,par.mod.crust.vpvsmax,20);
% No = histcounts2(posterior.vpvs,posterior.zmoh,Xk,Xh); No = No/maxgrid(No);
% Ni = histcounts2(prior.vpvs,prior.zmoh,Xk,Xh); Ni = Ni/maxgrid(Ni);
% contourf(midpts(Xk),midpts(Xh),No',[0:0.1:1],'linestyle','none');
% contour(midpts(Xk),midpts(Xh),Ni',[0.7:0.1:1],'--k','LineWidth',1.5);
plot(prior.vpvs,prior.zmoh,'.k','markersize',2)
plot(posterior.vpvs,posterior.zmoh,'.r','markersize',2)
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
     'xlim',[par.mod.crust.vpvsmin,par.mod.crust.vpvsmax],...
     'ylim',[par.mod.sed.hmin+par.mod.crust.hmin par.mod.sed.hmax+par.mod.crust.hmax])
title('H-K comparison','fontsize',16)


%% Crustal Vp/Vs
subplot(345), cla, hold on
X = midpts(linspace(par.mod.crust.vpvsmin,par.mod.crust.vpvsmax,20));
No = hist(posterior.vpvs,X)/posterior.Nstored;
Ni = hist(prior.vpvs,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'xlim',[par.mod.crust.vpvsmin,par.mod.crust.vpvsmax],'ylim',[0 axlim(gca,4)]) 
title('Crust Vp/Vs ratio','fontsize',16)

%% Crust radial anisotropy
if par.mod.crust.ximin~=par.mod.crust.ximax
subplot(346), cla, hold on
X = midpts(linspace(par.mod.crust.ximin,par.mod.crust.ximax,20));
No = hist(posterior.xicrust,X)/posterior.Nstored;
Ni = hist(prior.cxi,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'xlim',[par.mod.crust.ximin,par.mod.crust.ximax],'ylim',[0 axlim(gca,4)]) 
title('Crust Xi value','fontsize',16)
end

%% Mantle radial anisotropy
if par.mod.mantle.ximin~=par.mod.mantle.ximax
subplot(347), cla, hold on
X = midpts(linspace(par.mod.mantle.ximin,par.mod.mantle.ximax,20));
No = hist(posterior.ximant,X)/posterior.Nstored;
Ni = hist(prior.mxi,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'xlim',[par.mod.mantle.ximin,par.mod.mantle.ximax],'ylim',[0 axlim(gca,4)]) 
title('Mantle Xi value','fontsize',16)
end

%% Width of mantle NVG
% find depth indices with negative grad
nvgg = -0.004;
% posterior
Lnvg_post = zeros(posterior.Nstored,1);
for ii = 1:posterior.Nstored
nvind = find([diff(posterior.VSmantle(ii,:)')./diff(posterior.zatdep)<nvgg;0] & posterior.zatdep>posterior.zmoh(ii));
if isempty(nvind), Lnvg_post(ii)=nan; continue; end
a = diff(nvind);
b = find([a;inf]>1);
c = diff([0;b]);% length of sequences with nvgs
di1 = cumsum(c); % end points of sequences with nvgs
di0 = di1-c+1; % start points of sequences with nvgs
nvindm = [nvind(di0(c==max(c))):nvind(di1(c==max(c)))]';
Lnvg_post(ii) = diff(posterior.zatdep(nvindm([1,end])));
% figure(222);clf 
% plot(posterior.VSmantle(ii,:),posterior.zatdep,'b'); set(gca,'ydir','reverse'); hold on
% plot(posterior.VSmantle(ii,nvindm),posterior.zatdep(nvindm),'r'); 
end
% prior
Lnvg_pri = zeros(prior.Nstored,1);
for ii = 1:prior.Nstored
nvind = find([diff(prior.VSmantle(ii,:)')./diff(prior.zatdep)<nvgg;0] & prior.zatdep>prior.zmoh(ii));
if isempty(nvind), Lnvg_pri(ii)=nan; continue; end
a = diff(nvind);
b = find([a;inf]>1);
c = diff([0;b]);% length of sequences with nvgs
di1 = cumsum(c); % end points of sequences with nvgs
di0 = di1-c+1; % start points of sequences with nvgs
nvindm = [nvind(di0(c==max(c))):nvind(di1(c==max(c)))]';
Lnvg_pri(ii) = diff(prior.zatdep(nvindm([1,end])));
end

subplot(348), cla, hold on
X = midpts(linspace(0,200,40));
No = hist(Lnvg_post,X)/posterior.Nstored;
Ni = hist(Lnvg_pri,X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'xlim',[0,200],'ylim',[0 axlim(gca,4)]) 
title('max NVG width (km)','fontsize',16)


%% velocity at depth
X = midpts(linspace(par.mod.mantle.vsmin,par.mod.mantle.vsmax,20));
zdo = [55,67,90,100,120,150,180,250]; 
for iz = 1:length(zdo)
    izdo(iz,1) = crossing(prior.zatdep,[],zdo(iz));
    izdo(iz,2) = crossing(posterior.zatdep,[],zdo(iz));
end
for iz = 1:length(zdo)
subplot(3,8,16+iz), cla, hold on
No = hist(posterior.VSmantle(:,izdo(iz,2)),X)/posterior.Nstored;
Ni = hist(prior.VSmantle(:,izdo(iz,1)),X)/prior.Nstored;
bar(X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
% legend(num2str(model_summary.zmantle(:)),'location','northwest')
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
        'xlim',[par.mod.mantle.vsmin par.mod.mantle.vsmax],'ylim',[0 axlim(gca,4)],...
        'xtick',unique(round_level([par.mod.mantle.vsmin:0.1:par.mod.mantle.vsmax],0.3))) 
title(sprintf('Vs at %.0f km',prior.zatdep(izdo(iz))),'fontsize',16)
end


%% title
htit = title_custom([par.data.stadeets.sta,' ',par.data.stadeets.nwk],0.05,'fontweight','bold','fontsize',25);

if ifsave
    save2pdf(19,ofile,'/');
end


end


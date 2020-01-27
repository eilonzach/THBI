function plot_MODEL_SUMMARY( model_summary,par,ifsave,ofile)
% plot_MODEL_SUMMARY( model_summary,par,ifsave,ofile )
%   
% function to plot prior model distribution

if nargin<3 || isempty(ifsave)
    ifsave = 0; % default is not to save
end
if nargin<4 || isempty(ofile)
    ofile = 'figs/model_summary_fig.pdf';
end

try
    Nmods = model_summary.Npass;
catch
    Nmods = model_summary.Nstored;
end
    
    
    
%% plot prior
figure(90), clf, set(gcf,'pos',[43 141 1774 961]);
cls = get(groot,'defaultAxesColorOrder');

subplot(3,3,1)
[N,X] = hist(model_summary.zsed,20);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
set(gca,'fontsize',14), title('Sed thickness (km)','fontsize',16)

subplot(334)
[N,X] = hist(model_summary.zmoh,20);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
set(gca,'fontsize',14), title('Moho depth (km)','fontsize',16)

subplot(332),cla
[N,X] = hist([model_summary.VSsedtop,model_summary.VSsedbot],20);
hb = bar(X,N/Nmods,1,'facecolor','none'); 
for ii = 1:length(hb), hb(ii).EdgeColor = cls(ii,:);end
legend('top','bottom','location','north')
set(gca,'fontsize',14), title('Vs seds (km/s)','fontsize',16)

subplot(333),cla
[N,X] = hist([model_summary.VScrusttop,model_summary.VScrustbot],20);
hb = bar(X,N/Nmods,1,'facecolor','none'); 
for ii = 1:length(hb), hb(ii).EdgeColor = cls(ii,:);end
legend('top','bottom','location','north')
set(gca,'fontsize',14), title('Vs crust (km/s)','fontsize',16)


subplot(337)
[N,X] = hist(model_summary.VSmantle,20);
hb = bar(X,N/Nmods,1,'facecolor','flat'); xlim([3.5 4.9]);
legend(num2str(model_summary.zatdep(:)),'location','northwest')
set(gca,'fontsize',14), title('Vs at 50,100,150,...,300','fontsize',16)

subplot(335)
% Xh = linspace(par.mod.sed.hmin+par.mod.crust.hmin,par.mod.sed.hmax+par.mod.crust.hmax,30);
% Xk = linspace(par.mod.crust.vpvsmin,par.mod.crust.vpvsmax,20);
% N = histcounts2(model_summary.vpvs,model_summary.zmoh,Xk,Xh); 
% N = N/maxgrid(N);
% contourf(midpts(Xk),midpts(Xh),N',[0:0.1:1],'linestyle','none');
plot(model_summary.vpvs,model_summary.zmoh,'.k','markersize',1)
set(gca,'fontsize',14,'box','on','linewidth',1.5,'layer','top',...
     'xlim',[par.mod.crust.vpvsmin,par.mod.crust.vpvsmax],...
     'ylim',[par.mod.sed.hmin+par.mod.crust.hmin par.mod.sed.hmax+par.mod.crust.hmax])
title('H-K comparison','fontsize',16)

subplot(336)
[N,X] = hist(model_summary.fdVSmoh,20);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
set(gca,'fontsize',14), title('fractional dVs at moho (%)','fontsize',16)

if isfield(model_summary,'vpvs')
subplot(338)
[N,X] = hist(model_summary.vpvs,20);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
% xlim([1.7 1.9]);
set(gca,'fontsize',14), title('Crust Vp/Vs ratio','fontsize',16)
end

if isfield(model_summary,'ximant')
subplot(339)
[N,X] = hist(model_summary.ximant,0.9:0.01:1.1);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
xlim([0.95 1.05]);
set(gca,'fontsize',14), title('Mantle xi (radial anis.)','fontsize',16)
end


if ifsave
    save2pdf(90,ofile,'/');
end


end


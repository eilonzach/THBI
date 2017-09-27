function plot_MODEL_SUMMARY( model_summary,ifsave,ofile)
% plot_MODEL_SUMMARY( model_summary,ifsave,ofile )
%   
% function to plot prior model distribution

if nargin<2 || isempty(ifsave)
    ifsave = 0; % default is not to save
end
if nargin<3 || isempty(ofile)
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
legend(num2str(model_summary.zmantle(:)),'location','northwest')
set(gca,'fontsize',14), title('Vs at 50,100,150,...,300','fontsize',16)

subplot(335)
[N,X] = hist(model_summary.fdVSsed,20);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
set(gca,'fontsize',14), title('fractional dVs at sed/crust (%)','fontsize',16)

subplot(336)
[N,X] = hist(model_summary.fdVSmoh,20);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
set(gca,'fontsize',14), title('fractional dVs at moho (%)','fontsize',16)

if isfield(model_summary,'vpvs')
subplot(338)
[N,X] = hist(model_summary.vpvs,20);
hb = bar(X,N/Nmods,1,'facecolor','none','edgecolor','r');
% xlim([1.7 1.9]);
set(gca,'fontsize',14), title('Mantle Vp/Vs ratio','fontsize',16)
end

if ifsave
    save2pdf(90,ofile,'/');
end


end


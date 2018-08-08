function plot_FIG2_FIT_MODEL( final_model,posterior,prior,par,ifsave,ofile)
% plot_FIG2_FIT_MODEL( final_model,posterior,prior,par,ifsave,ofile)
%   


if nargin < 5 || isempty(ifsave)
    ifsave=false;
end

if nargin < 6 || isempty(ofile)
    ofile='fig2_FIT_MODEL';
end


%%  =========================  DATA FIT  =========================  
figure(62),clf,set(gcf,'pos',[015 76 950 529])
ax1 = axes('position',[0.07 0.09 .27 0.85]); hold on
ax2 = axes('position',[0.36 0.58 0.27 0.36]); hold on
ax3 = axes('position',[0.36 0.09 0.27 0.36]); hold on
ax4 = axes('position',[0.67 0.58 0.27 0.36]); hold on
ax5 = axes('position',[0.71 0.09 0.21 0.36]); hold on

vslims = [3.2 5.1];


%%  =========================  MODEL FIT  =========================  

%% Vs profile
Z = final_model.Z;

% plot(ax1,TRvs,TRZ,'-b','Linewidth',2.5);
fill(ax1,[final_model.VSsig1(:,1);...
          flipud(final_model.VSsig1(:,2))],...
          [Z;flipud(Z)],'-','Linewidth',1.5,'Facecolor',[0.7 0.7 0.7],'Edgecolor',[0.6 0.6 0.6]);
plot(ax1,final_model.VSav,Z,'-r','Linewidth',2);
plot(ax1,final_model.VSsig2,Z,'-','color',[0.4 0.4 0.4],'Linewidth',1);

Zmoh(1) = final_model.Zd(2).mu;
Zmoh(2) = final_model.Zd(2).std;
plot(ax1,vslims(1,:)-0.1,Zmoh(1)*[1 1],'--k','linewidth',2)
text(ax1,vslims(1,1)+0.1,Zmoh(1)+8,sprintf('$%.1f \\pm %.1f$',Zmoh),...
    'fontsize',17,'interpreter','latex')    

set(ax1,'ydir','reverse','fontsize',14,'ytick',[0:25:max(Z)],...
'color','none','xlim',vslims);
title(ax1,'\textbf{Vs model}','fontsize',20,'Interpreter','latex')
xlabel(ax1,'\textbf{Vs (km/s)}','fontsize',18,'Interpreter','latex')
ylabel(ax1,'\textbf{Depth (km)}','fontsize',18,'Interpreter','latex')

%% Prior/posterior histograms/coveriances

% 100 km Vs
ind=mindex(posterior.zatdep,100);
X = midpts(linspace(par.mod.mantle.vsmin,par.mod.mantle.vsmax,30));
No = hist(posterior.VSmantle(:,ind),X)/posterior.Nstored;
Ni = hist(prior.VSmantle(:,ind),X)/prior.Nstored;
bar(ax2,X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(ax2,X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(ax2,'fontsize',14), 
xlim(ax2,[par.mod.mantle.vsmin,par.mod.mantle.vsmax])
% title(ax2,'\textbf{ km}','fontsize',18,'interpreter','latex')
xlabel(ax2,sprintf('\\textbf{%.0f km Vs (km/s)}',posterior.zatdep(ind)),'fontsize',18,'Interpreter','latex')

% crustal VpVs
X = midpts(linspace(par.mod.crust.vpvsmin,par.mod.crust.vpvsmax,30));
No = hist(posterior.vpvs,X)/posterior.Nstored;
Ni = hist(prior.vpvs,X)/prior.Nstored;
bar(ax2,X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(ax2,X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(ax2,'fontsize',14), 
xlim(ax2,[par.mod.crust.vpvsmin,par.mod.crust.vpvsmax])
xlabel(ax2,'\textbf{Crustal $\mathbf{V_P/V_S}$ ratio}','fontsize',18,'Interpreter','latex')

par.sta = regexprep(par.sta,'_','\\_');
title(ax2,['\textbf{',par.sta,' ',par.nwk,'}'],'fontsize',22,'interpreter','latex')


% Moho depth
cla(ax3)
X = midpts(linspace(par.mod.crust.hmin+par.mod.sed.hmin,par.mod.crust.hmax+par.mod.sed.hmax,60));
No = hist(posterior.zmoh(:,end),X)/posterior.Nstored;
Ni = hist(prior.zmoh(:,end),X)/prior.Nstored;
bar(ax3,X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(ax3,X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(ax3,'fontsize',14), 
xlim(ax3,[par.mod.crust.hmin+par.mod.sed.hmin,par.mod.crust.hmax+par.mod.sed.hmax])
% xlim(ax3,[45 65])
% title(ax3,'\textbf{Moho depth}','fontsize',18,'interpreter','latex')
xlabel(ax3,'\textbf{Moho depth (km)}','fontsize',18,'Interpreter','latex')

% 125 km Vs
ind=mindex(posterior.zatdep,125);
cla(ax4)
X = midpts(linspace(par.mod.mantle.vsmin,par.mod.mantle.vsmax,30));
No = hist(posterior.VSmantle(:,ind),X)/posterior.Nstored;
Ni = hist(prior.VSmantle(:,ind),X)/prior.Nstored;
bar(ax4,X,No','facecolor',[0.9 0.1 0.1],'edgecolor','none','BarWidth',1);
bar(ax4,X,Ni','facecolor','none','edgecolor',[0.2 0.2 0.2],'BarWidth',1,'LineWidth',1.5);
set(ax4,'fontsize',14), 
xlim(ax4,[par.mod.mantle.vsmin,par.mod.mantle.vsmax])
% title(ax2,'\textbf{ km}','fontsize',18,'interpreter','latex')
xlabel(ax4,sprintf('\\textbf{%.0f km Vs (km/s)}',posterior.zatdep(ind)),'fontsize',18,'Interpreter','latex')

% Mantle Vs covariance
ind1 = mindex(posterior.zatdep,65);
ind2 = mindex(posterior.zatdep,82);
cla(ax5)
% plot(ax5,prior.VSmantle(:,ind1),prior.VSmantle(:,ind2),'o','linewidth',0.3,'color',0.3*[1 1 1])
plot(ax5,posterior.VSmantle(:,ind1),posterior.VSmantle(:,ind2),'or','markeredgecolor','none','markerfacecolor','r')
plot(ax5,[par.mod.mantle.vsmin,par.mod.mantle.vsmax],[par.mod.mantle.vsmin,par.mod.mantle.vsmax],'b--','linewidth',1.5)%, axis equal
set(ax5,'fontsize',14,'xlim',[3.9 4.65],'ylim',[3.9 4.65])%, axis equal
xlabel(ax5,sprintf('\\textbf{%.0f km Vs (km/s)}',posterior.zatdep(ind1)),'fontsize',18,'Interpreter','latex')
ylabel(ax5,sprintf('\\textbf{%.0f km Vs (km/s)}',posterior.zatdep(ind2)),'fontsize',18,'Interpreter','latex')

global TRUEmodel;
if isstruct(TRUEmodel)
    vs = TRUEmodel.VS;
    vp = TRUEmodel.VP;
    Z = TRUEmodel.z;
end
plot(ax1,vs,Z,':b','Linewidth',1.5);


% 
% % VsCB-VsMT covariance
% cla(ax5)
% plot(ax5,prior.VScrustbot,prior.VSmantle(:,1),'o','linewidth',0.3,'color',0.3*[1 1 1])
% plot(ax5,posterior.VScrustbot,posterior.VSmanttop,'or','markeredgecolor','none','markerfacecolor','r')
% set(ax5,'fontsize',14,'xlim',[3.2 4.65],'ylim',[3.2 4.65])%, axis equal
% xlabel(ax5,'\textbf{Lower crust Vs (km/s)}','fontsize',18,'Interpreter','latex')
% ylabel(ax5,'\textbf{Upper mantle Vs (km/s)}','fontsize',18,'Interpreter','latex')
% 
% % dVS_moh-VsMT covariance
% cla(ax5)
% plot(ax5,prior.fdVSmoh,prior.VSmantle(:,1),'o','linewidth',0.3,'color',0.3*[1 1 1])
% plot(ax5,posterior.fdVSmoh,posterior.VSmanttop,'or','markeredgecolor','none','markerfacecolor','r')
% set(ax5,'fontsize',14,'xlim',[0 23],'ylim',[3.8 4.8])%, axis equal
% xlabel(ax5,'\textbf{Moho dVs/Vs (\%)}','fontsize',18,'Interpreter','latex')
% ylabel(ax5,'\textbf{65 km Vs (km/s)}','fontsize',18,'Interpreter','latex')

% % dVS_moh-VsMT covariance
% cla(ax5)
% plot(ax5,prior.fdVSmoh,prior.vpvs,'o','linewidth',0.3,'color',0.3*[1 1 1])
% plot(ax5,posterior.fdVSmoh,posterior.vpvs,'or','markeredgecolor','none','markerfacecolor','r')
% set(ax5,'fontsize',14,'xlim',[8 23],'ylim',[1.7 1.9])%, axis equal
% xlabel(ax5,'\textbf{Lower crust Vs (km/s)}','fontsize',18,'Interpreter','latex')
% ylabel(ax5,'\textbf{Upper mantle Vs (km/s)}','fontsize',18,'Interpreter','latex')

% dVS_moh-VsMT covariance
% cla(ax5)
% plot(ax5,prior.fdVSmoh,prior.VSmantle(:,1),'o','linewidth',0.3,'color',0.3*[1 1 1])
% plot(ax5,posterior.fdVSmoh,posterior.VSmanttop,'or','markeredgecolor','none','markerfacecolor','r')
% set(ax5,'fontsize',14,'xlim',[0 23],'ylim',[3.8 4.8])%, axis equal
% xlabel(ax5,'\textbf{Moho dVs/Vs (\%)}','fontsize',18,'Interpreter','latex')
% ylabel(ax5,'\textbf{65 km Vs (km/s)}','fontsize',18,'Interpreter','latex')


pause(0.001)


if ifsave
    save2pdf(62,ofile,'/');
end

end


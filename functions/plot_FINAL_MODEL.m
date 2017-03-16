function plot_FINAL_MODEL( final_model,posterior,ifsave,ofile,ifinctru )
% plot_FINAL_MODEL( final_model,posterior,ifsave,ofile )
%
% function to plot final model (and how it compares to the true model)

if nargin<3 || isempty(ifsave)
    ifsave = 0; % default is not to save
end
if nargin<4 || isempty(ofile)
    ofile = [pwd,'/figs/final_model_fig.pdf'];
end
if nargin<5 || isempty(ifinctru)
    ifinctru=false;
end

Z = final_model.Z;
VSbest = final_model.VSbest;
VSstd = final_model.VSstd;
VPbest = final_model.VPbest;
VPstd = final_model.VPstd;
rhobest = final_model.rhobest;
rhostd = final_model.rhostd;

%% GET TARGET MODEL for comparison
if ifinctru
    global TRUEmodel
    TRvs = TRUEmodel.vs;
    TRvp = TRUEmodel.vp;
    TRrho = TRUEmodel.rho;
    TRZ = TRUEmodel.Z;
else
    TRZ = nan; TRvs=nan; TRvp=nan; TRrho=nan;
end


%% PLOT FINAL MODEL
figure(95); clf; set(gcf,'pos',[120 308 984 790])
ax1 = subplot(1,7,[1,2]); hold on
ax2 = subplot(1,7,[3,4]); hold on
ax3 = subplot(1,7,[5,6]); hold on
ax4 = subplot(1,7,7); hold on


plot(ax1,TRvs,TRZ,'-b','Linewidth',1.5);
plot(ax1,VSbest,Z,'-r','Linewidth',2);
plot(ax1,VSbest+2*VSstd,Z,'-','color',[0.4 0.4 0.4],'Linewidth',2);
plot(ax1,VSbest-2*VSstd,Z,'-','color',[0.4 0.4 0.4],'Linewidth',2);
set(ax1,'ydir','reverse','fontsize',14,'ytick',[0:25:max(Z)],...
    'color','none','xlim',[3.2 5.1]);
title(ax1,'Vs','fontsize',20)
xlabel(ax1,'Vs (km/s)','fontsize',16)
ylabel(ax1,'Depth (km)','fontsize',16)

plot(ax2,TRvp,TRZ,'-b','Linewidth',1.5);
plot(ax2,VPbest,Z,'-r','Linewidth',2);
plot(ax2,VPbest+2*VPstd,Z,'-','color',[0.4 0.4 0.4],'Linewidth',2);
plot(ax2,VPbest-2*VPstd,Z,'-','color',[0.4 0.4 0.4],'Linewidth',2);
set(ax2,'ydir','reverse','fontsize',14,'yticklabel','',...
    'ytick',[0:25:max(Z)],'color','none','xlim',[6 9.3]);
title(ax2,'Vp','fontsize',20)
xlabel(ax2,'Vp (km/s)','fontsize',16)
% ylabel(ax2,'Depth (km)','fontsize',16)

plot(ax3,TRrho,TRZ,'-b','Linewidth',1.5);
plot(ax3,rhobest,Z,'-r','Linewidth',2);
plot(ax3,rhobest+2*rhostd,Z,'-','color',[0.4 0.4 0.4],'Linewidth',2);
plot(ax3,rhobest-2*rhostd,Z,'-','color',[0.4 0.4 0.4],'Linewidth',2);
set(ax3,'ydir','reverse','fontsize',14,'yticklabel','',...
    'ytick',[0:25:max(Z)],'color','none','xlim',[2.5 3.9]);
title(ax3,'rho','fontsize',20)
xlabel(ax3,'rho (g/cc)','fontsize',16)
% ylabel(ax3,'Depth (km)','fontsize',16)

Zdh = midpts([0:0.5:max(Z)]); %Zdh = Zdh(1:3:end);
%moho
nm = hist(posterior.zmoh,Zdh);
plot(ax4,nm/sum(nm),Zdh,'k')
fill(ax4,nm/sum(nm),Zdh,[0.5 0.5 0.5])
% seds
if ~all(posterior.zsed==0)
ns = hist(posterior.zsed,Zdh);
plot(ax4,ns/sum(ns),Zdh,'k')
fill(ax4,ns/sum(ns),Zdh,[0.5 0.5 0.5])
end
set(ax4,'ydir','reverse','fontsize',14,'yticklabel','','ytick',[0:25:max(Z)],'color','none');
title(ax4,'Disc.','fontsize',20)
xlabel(ax4,'fraction','fontsize',16)

%% SAVE
if ifsave
    save2pdf(95,ofile,'/');
end

end


function plot_DATAFITS( trudata,savedata,par,ifsave,ofile)
%plot_DATAFITS( trudata,savedata,par,ifsave,ofile)
%   
% function to plot all accepted predicted data against the true data for
% the salient components

if nargin < 4 || isempty(ifsave)
    ifsave=false;
end

if nargin < 5 || isempty(ofile)
    ofile='true_vs_saved_data';
end

ps_xlims = [-10 30];
sp_xlims = [-30 10];


figure(60),clf,set(gcf,'pos',[015 576 915 529])
ax1 = axes('position',[0.07 0.57 0.43 0.38]);
ax2 = axes('position',[0.07 0.09 0.43 0.38]);
ax3 = axes('position',[0.58 0.10 0.40 0.82]);

%% Ps
if ~isempty(trudata.PsRF.ZRT)
axes(ax1), hold on, set(gca,'fontsize',13,'ytick',0)
plot(ax1,savedata.PsRF.tt(savedata.bestmods,:)',savedata.PsRF.R(savedata.bestmods,:)','color',[0.9 0.9 0.9],'linewidth',0.2), 
plot(ax1,trudata.PsRF.tt,trudata.PsRF.ZRT(:,2),'k','linewidth',2.5), 
xlabel(ax1,'\textbf{Time from P arrival} $\longrightarrow$','fontsize',18,'interpreter','latex')
ylabel(ax1,'\textbf{Radial disp.}','fontsize',18,'interpreter','latex')
xlim(ps_xlims), ylim(1.1*max(abs(trudata.PsRF.ZRT(:,2)))*[-1 1])
end

%% Sp
if ~isempty(trudata.SpRF.ZRT)
axes(ax2), hold on, set(gca,'fontsize',13,'xdir','reverse','ytick',0)
plot(ax2,savedata.SpRF.tt(savedata.bestmods,:)',savedata.SpRF.Z(savedata.bestmods,:)','color',[0.9 0.9 0.9],'linewidth',0.2), 
plot(ax2,trudata.SpRF.tt,trudata.SpRF.ZRT(:,1),'k','linewidth',2.5), 
xlabel(ax2,'$\longleftarrow$ \textbf{Time from S arrival} ','fontsize',18,'interpreter','latex')
ylabel(ax2,'\textbf{Vertical disp.}','fontsize',18,'interpreter','latex')
xlim(sp_xlims), ylim(1.1*max(abs(trudata.SpRF.ZRT(:,1)))*[-1 1])
end

%% SW
if ~isempty(trudata.SW.phV)
axes(ax3), hold on
plot(ax3,1./savedata.SW.periods(savedata.bestmods,:)',savedata.SW.phV(savedata.bestmods,:)','color',[0.9 0.9 0.9],'linewidth',0.2)
plot(ax3,1./trudata.SW.periods,trudata.SW.phV,'ok','linewidth',3,'markersize',10);
set(ax3,'fontsize',15,'xlim',[0 1./trudata.SW.periods(1)])
xlabel(ax3,'\textbf{Frequency (Hz)}','fontsize',18,'Interpreter','latex')
ylabel(ax3,'\textbf{Phase Velocity (km/s)}','fontsize',18,'Interpreter','latex')
else
delete(ax3) 
end


pause(0.001)


if ifsave
    save2pdf(60,ofile);
end

end


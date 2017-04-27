function axs = plot_TRUvsPRE_WAVEFORMS( trudata,predata,ifsave,ofile)
%plot_TRUvsPRE_WAVEFORMS( trudata,predata,ifsave,ofile )
%   
% function to plot predicted and true seismograms (Vertical and Radial)
% Assumes date in 3-column ZRT matrices with equal sample rate 

if nargin < 3 || isempty(ifsave)
    ifsave=false;
end

if nargin < 4 || isempty(ofile)
    ofile='true_vs_predicted_data';
end

xlims = [-10 30;-30 10]; %[P;S]


figure(58),clf,set(gcf,'pos',[015 576 1915 529])
ax1 = axes('position',[0.03 0.52 0.17 0.4]); hold on
ax2 = axes('position',[0.03 0.09 0.17 0.4]); hold on
ax3 = axes('position',[0.22 0.52 0.17 0.4]); hold on
ax4 = axes('position',[0.22 0.09 0.17 0.4]); hold on
ax5 = axes('position',[0.41 0.52 0.17 0.4]); hold on
ax6 = axes('position',[0.41 0.09 0.17 0.4]); hold on
ax7 = axes('position',[0.60 0.52 0.17 0.4]); hold on
ax8 = axes('position',[0.60 0.09 0.17 0.4]); hold on
ax9 = axes('position',[0.80 0.10 0.17 0.82]); hold on

axs=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9];

dtypes = fieldnames(predata);
dtypes(strcmp(dtypes,'SW')) = [];

%% SW
if ~isempty(predata.SW.phV)
    hp(1) = plot(ax9,trudata.SW.periods,trudata.SW.phV,'k.-','linewidth',3,'markersize',40);
    hp(2) = plot(ax9,predata.SW.periods,predata.SW.phV,'r.-','linewidth',1.5,'markersize',30);
    hl = legend(ax9,hp,'True','Pred','Location','SouthEast'); set(hl,'fontsize',15);
    set(ax9,'fontsize',15)
    xlabel(ax9,'Period (s)','fontsize',18)
    ylabel(ax9,'Phase Velocity (km/s)','fontsize',18)
    title(ax9, 'SW','fontsize',22)
else
    delete(ax9) 
end

%% RFs
for id = 1:length(dtypes)
    dtype = dtypes{id};
    xa1 = axs(mod(2*id+3,8)); % order [5,7,1,3]
    xa2 = axs(mod(2*id+3,8)+1); % order [6,8,2,4]
    if strcmp(dtype(1),'P'), ps=1;elseif strcmp(dtype(1),'S'), ps=2; end
    
    if ~isempty(predata.(dtype)) && ~isempty(predata.(dtype)(1).ZRT)
        for itr = 1:length(trudata.(dtype))
            plot(xa1,trudata.(dtype)(itr).tt,trudata.(dtype)(itr).ZRT(:,1),'k','linewidth',2.5)
            plot(xa1,predata.(dtype)(1).tt,predata.(dtype)(1).ZRT(:,1),'r','linewidth',1.5)
            plot(xa2,trudata.(dtype)(itr).tt,trudata.(dtype)(itr).ZRT(:,2),'k','linewidth',2.5)
            plot(xa2,predata.(dtype)(1).tt,predata.(dtype)(1).ZRT(:,2),'r','linewidth',1.5)
        end
        set(xa1,'xlim',xlims(ps,:),...
                'ylim',1.1*max(abs(trudata.(dtype)(1).ZRT(:,1)))*[-1 1],...
                'fontsize',13,'xticklabel',[])
        set(xa2,'xlim',xlims(ps,:),...
                'ylim',1.1*max(abs(trudata.(dtype)(1).ZRT(:,2)))*[-1 1],...
                'fontsize',13)
        title(xa1, regexprep(dtype,'_','-'),'fontsize',22)
        xlabel(xa2, sprintf('Time from %s arrival',dtype(1)),'fontsize',18)
    else
        delete(xa1),delete(xa2) 
    end
    
end

pause(0.001)


if ifsave
    save2pdf(58,ofile,'/');
end

% %% Ps
% if isfield(predata,'PsRF') && ~isempty(predata.PsRF(1).ZRT)
% axes(ax5), hold on, set(gca,'fontsize',13,'xticklabel',[])
% axes(ax6), hold on, set(gca,'fontsize',13)
% for itr = 1:length(predata.PsRF)
% plot(trudata.PsRF.tt,trudata.PsRF.ZRT(:,1),'k','linewidth',2.5), title('Ps','fontsize',22)
% plot(predata.PsRF.tt,predata.PsRF.ZRT(:,1),'r','linewidth',1.5), 
% plot(trudata.PsRF.tt,trudata.PsRF.ZRT(:,2),'k','linewidth',2.5), xlabel('Time from P arrival','fontsize',18)
% plot(predata.PsRF.tt,predata.PsRF.ZRT(:,2),'r','linewidth',1.5), 
% end
% xlim(ax5,ps_xlims), ylim(ax5,1.1*max(abs(trudata.PsRF.ZRT(:,1)))*[-1 1])
% xlim(ax6,ps_xlims), ylim(ax6,1.1*max(abs(trudata.PsRF.ZRT(:,2)))*[-1 1])
% else
% delete(ax5),delete(ax6) 
% end
% 
% %% Sp
% if isfield(predata,'SpRF') && ~isempty(predata.SpRF.ZRT)
% axes(ax7), hold on, set(gca,'fontsize',13,'xticklabel',[])
% plot(trudata.SpRF.tt,trudata.SpRF.ZRT(:,1),'k','linewidth',2.5), title('Sp','fontsize',22)
% plot(predata.SpRF.tt,predata.SpRF.ZRT(:,1),'r','linewidth',1.5),
% xlim(sp_xlims), ylim(1.1*max(abs(trudata.SpRF.ZRT(:,1)))*[-1 1])
% axes(ax8), hold on, set(gca,'fontsize',13)
% plot(trudata.SpRF.tt,trudata.SpRF.ZRT(:,2),'k','linewidth',2.5), xlabel('Time from S arrival','fontsize',18)
% plot(predata.SpRF.tt,predata.SpRF.ZRT(:,2),'r','linewidth',1.5), 
% xlim(sp_xlims), ylim(1.1*max(abs(trudata.SpRF.ZRT(:,2)))*[-1 1])
% else
% delete(ax7),delete(ax8) 
% end
% 
% 
% %% Ps_lo
% if isfield(predata,'PsRF_lo') && ~isempty(predata.PsRF_lo.ZRT)
% axes(ax1), hold on, set(gca,'fontsize',13,'xticklabel',[])
% plot(trudata.PsRF_lo.tt,trudata.PsRF_lo.ZRT(:,1),'k','linewidth',2.5), title('Ps-lo','fontsize',22)
% plot(predata.PsRF_lo.tt,predata.PsRF_lo.ZRT(:,1),'r','linewidth',1.5), 
% xlim(ps_xlims), ylim(1.1*max(abs(trudata.PsRF_lo.ZRT(:,1)))*[-1 1])
% axes(ax2), hold on, set(gca,'fontsize',13)
% plot(trudata.PsRF_lo.tt,trudata.PsRF_lo.ZRT(:,2),'k','linewidth',2.5), xlabel('Time from P arrival','fontsize',18)
% plot(predata.PsRF_lo.tt,predata.PsRF_lo.ZRT(:,2),'r','linewidth',1.5), 
% xlim(ps_xlims), ylim(1.1*max(abs(trudata.PsRF_lo.ZRT(:,2)))*[-1 1])
% else
% delete(ax1),delete(ax2) 
% end
% 
% %% Sp_lo
% if isfield(predata,'SpRF_lo') && ~isempty(predata.SpRF_lo.ZRT)
% axes(ax3), hold on, set(gca,'fontsize',13,'xticklabel',[])
% plot(trudata.SpRF_lo.tt,trudata.SpRF_lo.ZRT(:,1),'k','linewidth',2.5), title('Sp-lo','fontsize',22)
% plot(predata.SpRF_lo.tt,predata.SpRF_lo.ZRT(:,1),'r','linewidth',1.5),
% xlim(sp_xlims), ylim(1.1*max(abs(trudata.SpRF_lo.ZRT(:,1)))*[-1 1])
% axes(ax4), hold on, set(gca,'fontsize',13)
% plot(trudata.SpRF_lo.tt,trudata.SpRF_lo.ZRT(:,2),'k','linewidth',2.5), xlabel('Time from S arrival','fontsize',18)
% plot(predata.SpRF_lo.tt,predata.SpRF_lo.ZRT(:,2),'r','linewidth',1.5), 
% xlim(sp_xlims), ylim(1.1*max(abs(trudata.SpRF_lo.ZRT(:,2)))*[-1 1])
% else
% delete(ax3),delete(ax4) 
% end



end


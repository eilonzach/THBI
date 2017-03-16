function axs = plot_TRUvsPRE( trudata,predata,ifsave,ofile)
%plot_TRUvsPRE( trudata,predata,ifsave,ofile )
%   
% function to plot predicted and true seismograms (Vertical and Radial)
% Assumes date in 3-column ZRT matrices with equal sample rate 

if nargin < 3 || isempty(ifsave)
    ifsave=false;
end

if nargin < 4 || isempty(ofile)
    ofile='true_vs_predicted_data';
end

figure(57),clf,set(gcf,'pos',[015 576 1400 600])
ax1 = axes('position',[0.03 0.53 0.30 0.4]); hold on
ax2 = axes('position',[0.03 0.09 0.30 0.4]); hold on
ax3 = axes('position',[0.365 0.53 0.30 0.4]); hold on
ax4 = axes('position',[0.365 0.09 0.30 0.4]); hold on
ax5 = axes('position',[0.71 0.09 0.27 0.83]); hold on

axs=[ax1,ax2,ax3,ax4,ax5];

dtypes = fieldnames(predata);
dtypes(strcmp(dtypes,'SW')) = [];

%% RFs
for id = 1:length(dtypes)
    dtype = dtypes{id};
    xa = axs(mod(id+1,4)+1); % order [3,4,1,2]
    samprate = trudata.(dtype)(1).samprate;
    if ~isempty(predata.(dtype)(1).ZRT)
        for itr = 1:length(trudata.(dtype))
            itr2 = min([itr,length(predata.(dtype))]);
            cc_v1 = conv(trudata.(dtype)(itr).ZRT(:,1),predata.(dtype)(itr2).ZRT(:,2),'full'); % Vobs*Hpre
            cc_v2 = conv(trudata.(dtype)(itr).ZRT(:,2),predata.(dtype)(itr2).ZRT(:,1),'full'); % Hobs*Vpre
            cc_max = max(abs([cc_v1;cc_v2]));
            cc_t = [0:length(cc_v1)-1]'./samprate;
            plot(xa,cc_t,cc_v1,'k','linewidth',2.5),
            plot(xa,cc_t,cc_v2,'r','linewidth',1.5), 
        end
        if strcmp(dtype(1),'P')
            xlims = [0 1.5*trudata.(dtype)(1).nsamp./samprate];
        elseif strcmp(dtype(1),'S')
            xlims = length(cc_t)./samprate - [1.5*trudata.(dtype)(1).nsamp./samprate 0];
        end
        set(xa,'fontsize',13,'xlim',xlims,'ylim',1.1*cc_max*[-1 1])
        title(xa,regexprep(dtype,'_','-'),'fontsize',22,'pos',[mean(xlims),0.85*cc_max,0])
    else
        delete(xa) 
    end
end

% %% Ps
% if isfield(predata,'PsRF') && ~isempty(predata.PsRF.ZRT)
% Ps_v1 = conv(trudata.PsRF.ZRT(:,1),predata.PsRF.ZRT(:,2),'full'); % Vobs*Hpre
% Ps_v2 = conv(trudata.PsRF.ZRT(:,2),predata.PsRF.ZRT(:,1),'full'); % Hobs*Vpre
% Ps_max = max(abs([Ps_v1;Ps_v2]));
% Ps_t = [0:length(Ps_v1)-1]'./trudata.PsRF.samprate;
% plot(ax3,Ps_t,Ps_v1,'k','linewidth',2.5),
% plot(ax3,Ps_t,Ps_v2,'r','linewidth',1.5), 
% set(ax3,'fontsize',13,'xlim',ps_xlims,'ylim',1.1*Ps_max*[-1 1])
% title(ax3,'Ps','fontsize',22,'pos',[mean(ps_xlims),0.85*Ps_max,0])
% else
% delete(ax3) 
% end
% 
% %% Sp
% if isfield(predata,'SpRF') && ~isempty(predata.SpRF.ZRT)
% Sp_v1 = conv(trudata.SpRF.ZRT(:,1),predata.SpRF.ZRT(:,2),'full'); % Vobs*Hpre
% Sp_v2 = conv(trudata.SpRF.ZRT(:,2),predata.SpRF.ZRT(:,1),'full'); % Hobs*Vpre
% Sp_max = max(abs([Sp_v1;Sp_v2]));
% Sp_t = (length(Sp_v1)-[1:length(Sp_v1)]')./trudata.SpRF.samprate;
% plot(ax4,Sp_t,Sp_v1,'k','linewidth',2.5), 
% plot(ax4,Sp_t,Sp_v2,'r','linewidth',1.5), 
% set(ax4,'fontsize',13,'xlim',sp_xlims,'ylim',1.1*Sp_max*[-1 1])
% xlabel(ax4,'Time','fontsize',18)
% title(ax4,'Sp','fontsize',22,'pos',[mean(sp_xlims),0.85*Sp_max,0])
% else
% delete(ax4) 
% end
% 
% 
% %% Ps_lo
% if isfield(predata,'PsRF_lo') && ~isempty(predata.PsRF_lo.ZRT)
% Ps_lo_v1 = conv(trudata.PsRF_lo.ZRT(:,1),predata.PsRF_lo.ZRT(:,2),'full'); % Vobs*Hpre
% Ps_lo_v2 = conv(trudata.PsRF_lo.ZRT(:,2),predata.PsRF_lo.ZRT(:,1),'full'); % Hobs*Vpre
% Ps_lo_max = max(abs([Ps_lo_v1;Ps_lo_v2]));
% Ps_lo_t = [0:length(Ps_lo_v1)-1]'./trudata.PsRF.samprate;
% plot(ax1,Ps_lo_t,Ps_lo_v1,'k','linewidth',2.5)
% plot(ax1,Ps_lo_t,Ps_lo_v2,'r','linewidth',1.5), 
% set(ax1,'fontsize',13,'xlim',ps_xlims,'ylim',1.1*Ps_lo_max*[-1 1])
% ylabel(ax1,'Normalised Amp.','fontsize',18)
% title(ax1,'Ps-lo','fontsize',22,'pos',[mean(ps_xlims),0.85*Ps_lo_max,0])
% else
% delete(ax1) 
% ylabel(ax3,'Normalised Amp.','fontsize',18)
% end
% 
% %% Sp_lo
% if isfield(predata,'SpRF_lo') && ~isempty(predata.SpRF_lo.ZRT)
% Sp_lo_v1 = conv(trudata.SpRF_lo.ZRT(:,1),predata.SpRF_lo.ZRT(:,2),'full'); % Vobs*Hpre
% Sp_lo_v2 = conv(trudata.SpRF_lo.ZRT(:,2),predata.SpRF_lo.ZRT(:,1),'full'); % Hobs*Vpre
% Sp_lo_max = max(abs([Sp_lo_v1;Sp_lo_v2]));
% Sp_lo_t = (length(Sp_lo_v1)-[1:length(Sp_lo_v1)]')./trudata.SpRF_lo.samprate;
% plot(ax2,Sp_lo_t,Sp_lo_v1,'k','linewidth',2.5)
% plot(ax2,Sp_lo_t,Sp_lo_v2,'r','linewidth',1.5), 
% set(ax2,'fontsize',13,'xlim',sp_xlims,'ylim',1.1*Sp_lo_max*[-1 1])
% ylabel(ax2,'Normalised Amp.','fontsize',18)
% xlabel(ax2,'Time','fontsize',18)
% title(ax2,'Sp-lo','fontsize',22,'pos',[mean(sp_xlims),0.85*Sp_lo_max,0])
% else
% delete(ax2) 
% ylabel(ax4,'Normalised Amp.','fontsize',18)
% end

%% SW
if ~isempty(predata.SW.phV)
axes(ax5), hold on
hp(1) = plot(1./trudata.SW.periods,trudata.SW.phV,'k.-','linewidth',3,'markersize',40);
hp(2) = plot(1./predata.SW.periods,predata.SW.phV,'r.-','linewidth',1.5,'markersize',30);
hl = legend(hp,'True','Pred','Location','NorthEast'); set(hl,'fontsize',15);
set(ax5,'fontsize',15,'xlim',[0 1.1/min(trudata.SW.periods)])
xlabel('Frequency (Hz)','fontsize',18)
ylabel('Phase Velocity (km/s)','fontsize',18)
title(ax5,'SW','fontsize',22)
else
delete(ax5) 
end


pause(0.001)


if ifsave
    save2pdf(57,ofile,'/');
end

end


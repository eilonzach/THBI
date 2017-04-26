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
    if ~isempty(predata.(dtype)) && ~isempty(predata.(dtype).ZRT)
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


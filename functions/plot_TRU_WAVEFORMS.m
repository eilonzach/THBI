function axs = plot_TRU_WAVEFORMS( trudata,ifsave,ofile)
%plot_TRU_WAVEFORMS( trudata,ifsave,ofile )
%   
% function to plot predicted and true seismograms (Vertical and Radial)
% Assumes date in 3-column ZRT matrices with equal sample rate 

if nargin < 3 || isempty(ifsave)
    ifsave=false;
end

if nargin < 4 || isempty(ofile)
    ofile='true_data';
end

% xlims = [-5 26;-31 5]; %[P;S]



figure(58),clf,set(gcf,'pos',[15 54 1046 1051])
px = 0.28;
py = 0.12;
sy = 0.4;
axs = zeros(6,3);
for ix = 1:2
for iy = 1:3
for iyy = 1:2
    xl = 0.04 + 1.15*px*(ix-1);
    yc = 0.5 + 2.7*py*(2-iy);
    yb = -mod(iyy+1,2)*1.05*py + yc;
    axs(2*(iy-1)+iyy,ix) = axes('position',[xl yb px py]); hold on        
end
end
end
axs(1,3) = axes('position',[(0.05+2.3*px) (axpos(axs(1,1),2)+axpos(axs(1,1),4)-sy) px sy]); hold on
axs(2,3) = axes('position',[(0.05+2.3*px) axpos(axs(end,1),2) px sy]); hold on

% 
% ax1 = axes('position',[0.03 0.52 0.17 0.4]); hold on
% ax2 = axes('position',[0.03 0.09 0.17 0.4]); hold on
% ax3 = axes('position',[0.22 0.52 0.17 0.4]); hold on
% ax4 = axes('position',[0.22 0.09 0.17 0.4]); hold on
% ax5 = axes('position',[0.41 0.52 0.17 0.4]); hold on
% ax6 = axes('position',[0.41 0.09 0.17 0.4]); hold on
% ax7 = axes('position',[0.60 0.52 0.17 0.4]); hold on
% ax8 = axes('position',[0.60 0.09 0.17 0.4]); hold on
% ax9 = axes('position',[0.80 0.10 0.17 0.82]); hold on

% axs=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9];

dtypes = fieldnames(trudata);
axus = zeros(length(dtypes),2);
for id = 1:length(dtypes);
%     allpdytp(id,:)=parse_dtype(dtypes{id});
    pdtyp = parse_dtype(dtypes{id});
    if strcmp(pdtyp{1},'BW')
        if strcmp(pdtyp{2},'Ps')
            ix = 1;
        elseif strcmp(pdtyp{2},'Sp')
            ix = 2;
        end
        if strcmp(pdtyp{3},'def')
            iy = [1,2];
        elseif strcmp(pdtyp{3},'cms')
            iy = [5,6];
        end
        if strcmp(pdtyp{4},'lo')
            iy = [3,4];
        end
    elseif strcmp(pdtyp{1},'SW')
            ix = 3;
        if strcmp(pdtyp{2},'Ray')
            iy = 1;
        elseif strcmp(pdtyp{2},'Lov')
            iy = 2;
        end
    end
    axus(id,1:length(iy)) = axs(iy,ix);
%     title(axus(id,1),regexprep(dtypes{id},'_','-'))
end

delete(setdiff(axs,axus));


for id = 1:length(dtypes)
dtype = dtypes{id};
pdtyp = parse_dtype( dtype );
switch pdtyp{1}

%% SW
    case 'SW'
    
    hp(1) = plot(axus(id,1),trudata.(dtype).periods,trudata.(dtype).phV,'k.-','linewidth',3,'markersize',40);
    hl = legend(axus(id,1),hp,'True','Location','SouthEast'); set(hl,'fontsize',15);
    set(axus(id,1),'fontsize',15)
    xlabel(axus(id,1),'Period (s)','fontsize',18)
    ylabel(axus(id,1),'Phase Velocity (km/s)','fontsize',18)
    title(axus(id,1), 'SW','fontsize',22)

%% RFs
    case 'BW'

    xa1 = axus(id,1); % order [5,7,1,3]
    xa2 = axus(id,2); % order [6,8,2,4]
    if strcmp(pdtyp{2}(1),'P'), xp=1;xsv=0.2;elseif strcmp(pdtyp{2}(1),'S'), xp=0.2;xsv=1; end
    
    if ~isempty(trudata.(dtype)) && ~isempty(trudata.(dtype)(1).PSV)
        trunrm = zeros(length(trudata.(dtype)),1);
        for itr = 1:length(trudata.(dtype))
            trudata.(dtype)(itr).tt  = trudata.(dtype)(itr).tt(~isnan(trudata.(dtype)(itr).PSV(:,1)));
            trudata.(dtype)(itr).PSV = trudata.(dtype)(itr).PSV(~isnan(trudata.(dtype)(itr).PSV(:,1)),:);
%             trumax(itr) = max(abs(maxab(trudata.(dtype)(itr).PSV))); % get the max, to normalise trace
            trunrm(itr) = maxgrid(trudata.(dtype)(itr).PSV); % get the norm of the trace, to normalise the power
            plot(xa1,trudata.(dtype)(itr).tt,trudata.(dtype)(itr).PSV(:,1)./trunrm(itr),'k','linewidth',2.5)
            plot(xa2,trudata.(dtype)(itr).tt,trudata.(dtype)(itr).PSV(:,2)./trunrm(itr),'k','linewidth',2.5)
        end
        xlims = [min(trudata.(dtype)(itr).tt),max(trudata.(dtype)(itr).tt)];
        set(xa1,'xlim',xlims,...
                'ylim',1.1*xp*[-1 1],...
                'fontsize',13,'xticklabel',[])
        set(xa2,'xlim',xlims,...
                'ylim',1.1*xsv*[-1 1],...
                'fontsize',13)
        title(xa1, regexprep(dtype,'_','-'),'fontsize',22)
        xlabel(xa2, sprintf('Time from %s arrival',pdtyp{2}(1)),'fontsize',16)
    else
        delete(xa1),delete(xa2) 
    end

pause(0.001)


end % on switch

end % on dtype


if ifsave
    save2pdf(58,ofile,'/');
end


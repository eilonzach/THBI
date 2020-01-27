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

figure(57),clf,set(gcf,'pos',[015 576 1750 600])
ax1 = axes('position',[0.03 0.53 0.235 0.4]); hold on
ax2 = axes('position',[0.29 0.53 0.24 0.4]); hold on
ax3 = axes('position',[0.03 0.09 0.235 0.4]); hold on
ax4 = axes('position',[0.29 0.09 0.24 0.4]); hold on
ax5 = axes('position',[0.565 0.09 0.205 0.83]); hold on
ax6 = axes('position',[0.795 0.09 0.19 0.83]); hold on

% ax11 = axes('position',[0.03 0.09 0.11 0.84]); hold on
% ax12 = axes('position',[0.155 0.09 0.11 0.84]); hold on

if any(strcmp(fieldnames(trudata),'HKstack_P'))
    ax7 = axes('pos',[axpos(ax3,[1,2,3]) sum(axpos(ax1,[2,4]))-axpos(ax3,2)]); hold on
    delete([ax1,ax3]);
end


axs=[ax1,ax2,ax3,ax4,ax5,ax6];

dtypes = fieldnames(predata);

for id = 1:length(dtypes)
dtype = dtypes{id};
[ pdtyp ] = parse_dtype( dtype );
switch pdtyp{1}

%% BWs
    case {'BW','RF'}

    if strcmp(pdtyp{2},'Ps'), ipl = 3; elseif strcmp(pdtyp{2},'Sp'), ipl = 4; end
    if ~strcmp(pdtyp{3},'def') || ~strcmp(pdtyp{4},'def'), ipl = ipl-2; end
    xa = axs(ipl); % order [3,4,1,2]
    
    if strcmp(pdtyp{3},'ccp'), trudata.(dtype)(1).samprate = 1;end
    
    samprate = trudata.(dtype)(1).samprate;
    if ~isempty(predata.(dtype)) && ~isempty(predata.(dtype)(1).PSV)
        for itr = 1:length(trudata.(dtype))
%             itr2 = min([itr,length(predata.(dtype))]);
            cc_v1 = conv(trudata.(dtype)(itr).PSV(:,1),predata.(dtype)(itr).PSV(:,2),'full'); % Vobs*Hpre
            cc_v2 = conv(trudata.(dtype)(itr).PSV(:,2),predata.(dtype)(itr).PSV(:,1),'full'); % Hobs*Vpre
            cc_max = max(abs([cc_v1;cc_v2]));
            cc_t = [0:length(cc_v1)-1]'./samprate;
            plot(xa,cc_t,cc_v1,'k','linewidth',2.5),
            plot(xa,cc_t,cc_v2,'r','linewidth',1.5), 
        end
        if strcmp(pdtyp{2}(1),'P')
            xlims = [0 1.5*trudata.(dtype)(1).nsamp./samprate];
        elseif strcmp(pdtyp{2}(1),'S') && ~strcmp(pdtyp{3},'ccp')
            xlims = length(cc_t)./samprate - [1.5*trudata.(dtype)(1).nsamp./samprate 0];
        elseif strcmp(pdtyp{2}(1),'S') && strcmp(pdtyp{3},'ccp')
            xlims = [0 600];
        end
        set(xa,'fontsize',13,'xlim',xlims,'ylim',max([1.1*cc_max,axlim(xa,4)])*[-1 1])
        title(xa,regexprep(dtype,'_','-'),'fontsize',22,'pos',[mean(xlims),0.85*cc_max,0])
    else
        delete(xa) 
    end


%% SW
    case 'SW'
        switch pdtyp{2}
            
            case {'Ray','Lov'}
            
            if strcmp(pdtyp{2},'Ray'),msym = 'x'; else msym = 'o'; end
                
            axes(ax5), hold on
            hp(1) = plot(1./trudata.(dtype).periods,trudata.(dtype).phV,'k-','linewidth',3,'markersize',12,'Marker',msym);
            hp(2) = plot(1./predata.(dtype).periods,predata.(dtype).phV,'r-','linewidth',1.5,'markersize',10,'Marker',msym);
            hl = legend(hp,'True','Pred','Location','NorthEast'); set(hl,'fontsize',15);
            figN_add('x=R, o=L',ax5,0.05,0.04,12);
            set(ax5,'fontsize',15,'xlim',[0 1.1/min(trudata.(dtype).periods)])
            xlabel('Frequency (Hz)','fontsize',18)
            ylabel('Phase Velocity (km/s)','fontsize',18)
            title(ax5,'SW','fontsize',22)
            

            case 'HV' 
    
            axes(ax6), hold on
            errorbar(1./trudata.(dtype).periods,trudata.(dtype).HVr,2*trudata.(dtype).sigma.*ones(size(trudata.(dtype).periods)),'k')
            hp(1) = plot(1./trudata.(dtype).periods,trudata.(dtype).HVr,'k.-','linewidth',3,'markersize',40);
            hp(2) = plot(1./predata.(dtype).periods,predata.(dtype).HVr,'r.-','linewidth',1.5,'markersize',30);
            hl = legend(hp,'True','Pred','Location','NorthEast'); set(hl,'fontsize',15);
            set(ax6,'fontsize',15,'xlim',[0 1.1/min(trudata.(dtype).periods)])
            xlabel('Frequency (Hz)','fontsize',18)
%             ylabel('HV ratio)','fontsize',18)
            title(ax6,'HV ratio','fontsize',22)
        end
   
%% HKstack        
    case {'HKstack'}
        contourf(ax7,trudata.(dtype).K,trudata.(dtype).H,trudata.(dtype).Esum',30,'linestyle','none');
        colorbar(ax7,'southoutside')
        
        plot(predata.HKstack_P.K,predata.HKstack_P.H,'ok','linewidth',2,...
            'markerfacecolor','r','markersize',7)
        
        title(ax7, regexprep(dtype,'_','-'),'fontsize',22)
        xlabel(ax7, 'Vp/Vs ratio','fontsize',16)
        ylabel(ax7, 'Moho depth','fontsize',16)
        set(ax7,'fontsize',13,'ydir','reverse')
        
        

end


pause(0.001)


if ifsave
    save2pdf(57,ofile,'/');
end

end


function plot_KNOT_TRENDS( allmodels,par,ofile )
%plot_KNOT_TRENDS( allmodels,par,ofile )

if nargin<3
    ofile=[];
end

downsampfac = 1;

figure(83), clf, set(gcf,'pos',[65 297 1127 808])

ax1 = axes('pos',[0.07,0.2,0.73,0.7]); hold on
ax2 = axes('pos',[0.85,0.2,0.1,0.7]); hold on
ax3 = axes('pos',[0.07,0.07,0.73,0.1]); hold on
ax4 = axes('pos',[0.85,0.07,0.1,0.1]); hold on

%% loop through each chain
if iscell(allmodels),nchains = length(allmodels); else nchains=1; end
for iii = 1:nchains
    basecol = colour_get(iii,nchains+1,0,parula); basecol = basecol(:)';
    
    if nchains>1 
        am = allmodels{iii};
        if isempty(am), continue; end
    else 
        am = allmodels;
        if iscell(am), am = am{1}; end
    end
    
    downsamp = (1:downsampfac:length(am));
    ydownsamp = false(length(am),1); ydownsamp(downsamp)=true;

    %% Knot locations
    knmat_c = nan(par.mod.crust.kmax,length(am));
    knmat_m = nan(par.mod.mantle.kmax,length(am));

    for im = 1:length(am)
        knmat_c(1:am(im).crustmparm.Nkn,im) = am(im).crustmparm.knots;
        knmat_m(1:am(im).mantmparm.Nkn,im) = am(im).mantmparm.knots;
    end
    
    itmat_c = ones(par.mod.crust.kmax,1)*[am.iter];
    itmat_m = ones(par.mod.mantle.kmax,1)*[am.iter];
    
    %% Histogram of node locations
    knmat_c_use = knmat_c(2:end,[am.bestmods]);
    knmat_m_use = knmat_m(2:end,[am.bestmods]);
    knzc = reshape(knmat_c_use,numel(knmat_c_use),1);
    knzm = reshape(knmat_m_use,numel(knmat_m_use),1);
    knzm(isnan(knzm))=[];knzm(knzm==par.mod.maxz)=[];
    knzc(isnan(knzc))=[];knzc(knzc==0)=[];
    knz = [knzc;knzm];
    Xkn = midpts([0:par.mod.dz:par.mod.maxz]);
    pkn = hist(knz,Xkn)/sum([am.bestmods]);
    % fill out edges w/ zeros
    Xkn = [0,Xkn,par.mod.maxz]'; pkn = [0,pkn,0]';
    
    %% # of splines
    Nsp = zeros(length(am),1);
    for ii = 1:length(am)
        Nsp(ii) = am(ii).crustmparm.Nsp + am(ii).mantmparm.Nsp;
    end
    Xnsp = midpts([0:1:par.mod.mantle.kmax+par.mod.crust.kmax+2]);
    pnsp = hist(Nsp([am.bestmods]),Xnsp)/sum([am.bestmods]);
    knlim = 2+[par.mod.mantle.kmin+par.mod.crust.kmin,par.mod.mantle.kmax+par.mod.crust.kmax]; if diff(knlim)==0, knlim = knlim+[-0.5 0.5]; end
    
    %% plotting

    % knot locations
    plot(ax1,itmat_c(:,downsamp),knmat_c(:,downsamp),'.','color',[0.8 0.8 0.8])
    plot(ax1,itmat_c(:,[am.bestmods]'&ydownsamp),knmat_c(:,[am.bestmods]'&ydownsamp),'s','color',basecol,'markerfacecolor',basecol)
    plot(ax1,itmat_m(:,downsamp),knmat_m(:,downsamp),'.','color',[0.8 0.8 0.8])
    plot(ax1,itmat_m(:,[am.bestmods]'&ydownsamp),knmat_m(:,[am.bestmods]'&ydownsamp),'o','color',basecol,'markerfacecolor',basecol)
    plot(ax1,[am.iter],[am.zmoh],'-k','linewidth',2)
    plot(ax1,[am.iter],[am.zsed],'--k','linewidth',2)
    set(ax1,'ydir','reverse','fontsize',15,'xticklabel',[])
    title(ax1,'Knot placement','fontsize',22)
    ylabel(ax1,'Depth (km)','fontsize',20)
   
    % pdf of knot locations

    fill(ax2,pkn,Xkn,basecol)
    plot(ax2,pkn,Xkn,'k')
    set(ax2,'ydir','reverse','fontsize',15,'xlim',[0,1.05*max(pkn)])
    
    % # of splines
    plot(ax3,[am.iter],Nsp,'.','color',[0.8 0.8 0.8])
    plot(ax3,[am([am.bestmods]).iter],Nsp([am.bestmods]),'o','color',basecol,'markerfacecolor',basecol)
    set(ax3,'fontsize',15,'ylim',knlim) 
    xlabel(ax3,'Iteration','fontsize',20)
    
    % pdf of splines
    fill(ax4,pnsp,Xnsp,basecol)
    plot(ax4,pnsp,Xnsp,'k')
    set(ax4,'fontsize',15,...
        'ylim',knlim,...
        'xlim',[0,1.05*max(pnsp)]) 
    
end 


if ~isempty(ofile)
    if nchains>1
        save2jpg(83,ofile,'/')
    else
        save2pdf(83,ofile,'/')
    end
end

end




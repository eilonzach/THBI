function plot_invtime(misfits_perchain,ofile)
%plot_invtime(misfits,ofile)
%   Plot the time taken by the inversion

if nargin<2
    ofile=[];
end

downsampfac = 4;

figure(66), clf, set(gcf,'pos',[1654 719 883 618]), 
ax1 = axes('pos',[0.07 0.54 0.88 0.4]); hold on
ax2 = axes('pos',[0.07 0.08 0.88 0.4]); hold on



%% loop through each chain
if iscell(misfits_perchain),nchains = length(misfits_perchain); else nchains=1; end
for iii = 1:nchains
    basecol = colour_get(iii,nchains+1,0,parula); basecol = basecol(:)';
    
    if nchains>1
        mf = misfits_perchain{iii};
        if isempty(mf); continue; end
    else 
        mf = misfits_perchain; 
        if iscell(mf), mf = mf{1}; end
    end
    
    downsamp = (1:downsampfac:length(mf.iter));
    
    iters = mf.iter; iters(length(mf.iter)+1:length(mf.time)) = nan; 
    iters = iters(downsamp);
    cumtimes = mf.time(downsamp);
    avgtimes = diff([0;cumtimes])./diff([0;iters]);
    
    % get time unit for cumtimes
    timestr = {'s','min','hr','day'}; 
    timenrm = [1,60,60*60,24*60*60];
    for it = 1:length(timestr)
        if max(cumtimes)>timenrm(it), ituse = it; end
    end
        

    plot(ax1,iters,cumtimes./timenrm(ituse),'-o','linewidth',1.5,'color',basecol);
    plot(ax2,iters,avgtimes,'-o','linewidth',1.5,'color',basecol);
    plot(ax2,[0 max(mf.iter)],(max(mf.time)./max(mf.iter))*[1 1],'--r','linewidth',1.5);
    
    set(ax1,'xlim',[0 max(iters)],'fontsize',16)
	set(ax2,'xlim',[0 max(iters)],'fontsize',16)
end
    
xlabel(ax2,'Iteration','fontsize',20,'interpreter','latex')


ylabel(ax1,sprintf('\\textbf{Cumulative time (%s)}',timestr{ituse}),'fontsize',20,'interpreter','latex')
ylabel(ax2,'\textbf{Avg time per iter (s)}','fontsize',20,'interpreter','latex')

title(ax1,'\textbf{Time taken by inversion}','fontsize',22,'interpreter','latex')     

if ~isempty(ofile)
    save2pdf(66,ofile,'/')
end


end


function plot_clust_radplot( seazs,gcarcs,iclust,ax,ifsave,ofile )
%plot_clust_radplot( seazs,gcarcs,iclust )
%   Plot the events by cluster on a radial plot

if nargin<4 || isempty(ax)
    figure(34), clf, set(gcf,'pos',[202 53 700 700])
    ax = polaraxes('ThetaZeroLocation','top');
end
if nargin<5 || isempty(ifsave)
    ifsave = false;
end
if nargin<6 || isempty(ofile)
    ofile = './evt_dibn.pdf';
end

minarc = 25;
maxarc = 75;

nclust = length(unique(iclust));


for icl = 1:nclust
try
    hp = polarscatter(ax,d2r(seazs(iclust==icl)),(gcarcs(iclust==icl)-minarc),30,colour_get(icl,nclust,1,parula)','filled','markeredgecolor','k');
end
    hold on
end
plot(ax,0,0,'^k','linewidth',2,'markersize',13,'markerfacecolor',[0.6 0.6 0.6])
    
set(ax,'rtick',[10:10:maxarc-minarc],'rticklabel',[minarc+10:10:maxarc])
set(ax,'linewidth',2,'fontsize',17,'ThetaZeroLocation','top','ThetaDir','clockwise')

if ifsave
    save2pdf(gcf,ofile)
end


end


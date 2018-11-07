function [ iclust,nclust ] = eqcluster( seazs,gcarcs,min4clust,bazwin,arcwin,ifplot )
%EQCLUSTER Summary of this function goes here
%   Detailed explanation goes here

if nargin<3 || isempty(min4clust)
    min4clust = 10;
end
if nargin<4 || isempty(bazwin)
    bazwin = 10;
end
if nargin<5 || isempty(arcwin)
    arcwin = 15;
end
if nargin<6 || isempty(ifplot)
    ifplot = false;
end

iclust = zeros(size(seazs));

%% Parse EQ into gcarc/baz space
baztry = [0:359]';
arctry = [00:180-arcwin+1]';
eqinbin = zeros(length(baztry),length(arctry));

for ibaz = 1:length(baztry)
    for iarc = 1:length(arctry)
        baz = baztry(ibaz);
        arc = arctry(iarc);
        eqinbin(ibaz,iarc) = sum(abs(azdiff(baz,seazs))<bazwin/2 & abs(gcarcs-arc)<arcwin/2);
    end
end

eqinbin_full = eqinbin; 
eqinbin_full(eqinbin_full==0) = nan;

%% recursively identify clusters
icl = 0;
while any(any(eqinbin>min4clust))
    icl = icl+1;
    [~,iarc,ibaz] = maxgrid(eqinbin);
    bazcl(icl) = baztry(ibaz);
    arccl(icl) = arctry(iarc);
    indclust = (abs(bazcl(icl)-seazs) <= bazwin/2) & (abs(arccl(icl)-gcarcs) <= arcwin/2);
    iclust(indclust) = icl;
    
    for ibaz = 1:length(baztry)
    for iarc = 1:length(arctry)
        if (abs(baztry(ibaz) - bazcl(icl)) < bazwin) & (abs(arctry(iarc) - arccl(icl)) <= arcwin)
            eqinbin(ibaz,iarc) = 0;

        end
    end
    end 
end
nclust = icl;

if nclust == 0 
    fprintf('Not enough earthquakes in any bin to form a single cluster\n');
    return
end

%% plot histograms of parms within each cluster
if ifplot
bazbins = [0:2:360]';
arcbins = [0:2:180]';
bazhist = zeros(length(bazbins),nclust);
archist = zeros(length(arcbins),nclust);
for icl = 1:nclust
    bazhist(:,icl) = hist(seazs(iclust==icl) ,bazbins);
    archist(:,icl) = hist(gcarcs(iclust==icl),arcbins);
end


figure(33), clf, set(gcf,'pos',[202 53 940 640])
ax1 = axes('pos',[0.3 0.3 0.65 0.65]);
ax2 = axes('pos',[0.08 0.3 0.2 0.65]);
ax3 = axes('pos',[0.3 0.08 0.65 0.2]);

contourf(ax1,baztry,arctry,eqinbin_full',100,'edgecolor','none')
bar(ax3,bazbins,bazhist,1,'stack'), colormap(ax2,jet)
barh(ax2,arcbins,archist,1,'stack'), colormap(ax3,jet)

set(ax1,'xlim',[0,360],'ylim',[20,max(gcarcs)+2],'xtick',[0:45:360],'xticklabel',[],'ytick',[0:10:180],'yticklabel',[])
set(ax3,'xlim',[0,360],'xtick',[0:45:360])
set(ax2,'ylim',[20 max(gcarcs)+2],'ytick',[0:10:180])

ylabel(ax2,'\textbf{Station-event distance} ($^{\circ}$)','interpreter','latex','fontsize',22)
xlabel(ax3,'\textbf{Station-event azimuth} ($^{\circ}$)','interpreter','latex','fontsize',22)


end

end % on function
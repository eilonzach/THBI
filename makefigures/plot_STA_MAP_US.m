function [figh,axh,hstas] = plot_STA_MAP_US(fig,stainfo,lolim,lalim)

if nargin<4
    lalim = [39 51.5];
end 
if nargin<3
    lolim = [-130 -85];
end

addpath('~/Dropbox/MATLAB/lib/m_map/');

figh = figure(fig); clf; set(gcf,'pos',[70 708 1092 598]);

m_proj('Gall-peters','lon',lolim,'lat',lalim)
m_coast;
m_grid('linestyle','none','tickdir','out','linewidth',3,'fontsize',14,'layer','top');



% colormap(flipud(copper));
[CS,CH]=m_etopo2('contourf',[[-5000:200:0],[1,10,50,100:200:3000]],'edgecolor','none');
colormap(1 - (1-[m_colmap('water',27);m_colmap('bland',18)]).^1);
caxis([-5000 3000])
% [Cw,Cw]=m_etopo2('contourf',[0 0],'edgecolor','none');

% US state borders
addpath('~/Dropbox/MATLAB/lib/borders_v3.1.2/borders/');
[bordlat,bordlon] = borders('Continental US');
for ib = 1:length(bordlat)
    hbord(ib) = m_plot(bordlon{ib},bordlat{ib},'k');
end
set(hbord,'linewidth',1.3);

stinbounds = stainfo.slats<=max(lalim) & stainfo.slats>=min(lalim) &...
             stainfo.slons<=max(lolim) & stainfo.slons>=min(lolim);


hstas = m_plot(stainfo.slons(stinbounds),stainfo.slats(stinbounds),'^k','MarkerFaceColor','r','markersize',11,'linewidth',2);

axh = gca;



end

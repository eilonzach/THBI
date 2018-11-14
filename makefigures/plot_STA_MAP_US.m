% function figh = plot_STA_MAP_US(stainfo)

addpath('~/Dropbox/MATLAB/lib/m_map/');

figure(1); clf; set(gcf,'pos',[70 708 1092 598]);

m_proj('Gall-peters','lon',[-130 -85],'lat',[39 56])
m_coast;
m_grid('linestyle','none','tickdir','out','linewidth',3);
% m_elev('contourf',[0:400:6000]);
% colormap(flipud(copper));
[CS,CH]=m_etopo2('contourf',[-5000:200:3000],'edgecolor','none');
colormap(1 - (1-[m_colmap('blues',25);m_colmap('gland',15)]).^1);

m_plot(stainfo.slons,stainfo.slats,'^k','MarkerFaceColor','r','markersize',11,'linewidth',2)   


% end

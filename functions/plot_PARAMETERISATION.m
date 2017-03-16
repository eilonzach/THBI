function plot_PARAMETERISATION( model,ofile,ax )
% plot_PARAMETERISATION( model,ofile,ax )
% 
%  Function to plot the parameterisation of the model, including the
%  splines. If you give an ofile it will save the figure to that file.
%  Also, option not to make a new figure but to plot into custom axes (ax).
if nargin < 2|| isempty(ofile), ifsave = false; else ifsave = true; end
if nargin < 3 || isempty(ax), ax = []; end

z_sed = model.z(1:2);
z_crust = model.z(2+[1:size(model.crustmparm.splines,1)]);
z_mantle = model.z(length(z_crust)+1+2:end);

%% Make fig
if isempty(ax)
figure(3),clf,set(gcf,'pos',[100 100 400 800])
ax = gca;
end
cla(ax);
hold(ax,'on')

% SPLINES
% sed nodes
plot(ax,[1,0],z_sed,'k','markersize',15,'linewidth',3)
plot(ax,[0,1],z_sed,'k','markersize',15,'linewidth',3)
% crust splines
for is = 1:size(model.crustmparm.splines,2)
    plot(ax,model.crustmparm.splines(:,is),z_crust,'-g','linewidth',3)
end
% mantle splines
for is = 1:size(model.mantmparm.splines,2)
    plot(ax,model.mantmparm.splines(:,is),z_mantle,'-b','linewidth',3)
end

% BOUNDARIES
% surface
plot(ax,[0,1],max(z_sed)*[1,1],':','linewidth',3,'color',[0.5 0.5 0.5])
% sed/crust boundary
plot(ax,[0,1],max(z_sed)*[1,1],':','linewidth',3,'color',[0.5 0.5 0.5])
text(ax,0.5,max(z_sed),'Sed/crust','HorizontalAlignment','center','verticalalignment','top','fontsize',15,'fontweight','bold','color',[0.5 0.5 0.5])
% crust/mantle boundary
plot(ax,[0,1],max(z_crust)*[1,1],':','linewidth',3,'color',[0.5 0.5 0.5])
text(ax,0.5,max(z_crust),'Moho','HorizontalAlignment','center','verticalalignment','top','fontsize',15,'fontweight','bold','color',[0.5 0.5 0.5])

set(ax,'ydir','reverse','fontsize',14,'box','on','layer','top','linewidth',3,'xtick',[0 1])
title(ax,'Model parameterisation','fontsize',20)
xlabel(ax,'Weight','fontsize',17,'fontweight','bold')
ylabel(ax,'Depth (km)','fontsize',17,'fontweight','bold')

if ifsave
    save2pdf(gcf,ofile);
end


end


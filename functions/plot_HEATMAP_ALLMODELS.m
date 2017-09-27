function [hmps,hmpp] = plot_HEATMAP_ALLMODELS(suite_of_models,par,ifsave,ofile)
% [hmps,hmpp] = plot_HEATMAP_ALLMODELS(suite_of_models,par,ifsave,ofile)
% 
% function to plot a heatmap pdf of the velocity models

if nargin<3 || isempty(ifsave)
    ifsave = false; % default is not to save
end
if nargin<4 || isempty(ofile)
    ofile = 'figs/heatmap_models_fig.pdf';
end

%% GET TARGET MODEL for comparison
% global TRUEmodel
% TRvs = TRUEmodel.vs;
% TRvp = TRUEmodel.vp;
% TRrho = TRUEmodel.rho;
% TRZ = TRUEmodel.Z;

sm = suite_of_models;

Xvs = [par.mod.sed.vsmin:0.04:par.mod.mantle.vsmax];
Xvp = [1.3*par.mod.sed.vsmin:0.04:par.mod.mantle.vsmax*2];
Nz  = length(sm(1).Z);

hmps = zeros(Nz,length(Xvs));
hmpp = zeros(Nz,length(Xvp));
nsaved = 0;

% loop over depths, adding 
for iz = 1:Nz
    hmps(iz,:) = hist(sm.VS(iz,:),Xvs);
    hmpp(iz,:) = hist(sm.VP(iz,:),Xvp);
end
% norm pdf to 1
hmps = hmps/size(sm.VS,2);
hmpp = hmpp/size(sm.VS,2);
% log pdf
hmps = log(hmps); hmps(isinf(hmps)) = -20;
hmpp = log(hmpp); hmpp(isinf(hmpp)) = -20;


%% Plot
figure(25); clf; set(gcf,'pos',[120 151 920 947])
ax1 = axes(gcf,'pos',[0.13 0.11 0.35 0.815]); hold on
ax2 = axes(gcf,'pos',[0.52 0.11 0.35 0.815]); hold on
ax3 = axes(gcf,'pos',[0.91 0.11 0.02 0.815]); hold on


%% data
contourf(ax1,Xvs,sm.Z,hmps,[-5:0.1:-0.1],'edgecolor','none');
contourf(ax2,Xvp,sm.Z,hmpp,[-5:0.1:-0.1],'edgecolor','none');

%% pretty
xlim(ax1,[3 4.9])
xlim(ax2,[5.4 8.8])

set(ax1,'ydir','reverse','fontsize',15,'ytick',[0:25:max(sm.Z)],'color','none');
set(ax2,'ydir','reverse','fontsize',15,'yticklabel','','ytick',[0:25:max(sm.Z)],'color','none');
set(ax3,'visible','off')

title(ax1,'Vs','fontsize',24)
title(ax2,'Vp','fontsize',24)

xlabel(ax1,'Vs (km/s)','fontsize',18)
ylabel(ax1,'Depth (km)','fontsize',18)

xlabel(ax2,'Vp (km/s)','fontsize',18)

% colours
caxis(ax1,[-5 0])
caxis(ax2,[-5 0])
caxis(ax3,[-5 0])
hcb = colorbar(ax3);
set(hcb,'pos',[0.89 0.71 0.02 0.215],'fontsize',14,'yaxislocation','right')
hcby = ylabel(hcb,'log_{10}(Probability)','fontweight','bold','verticalalignment','middle','HorizontalAlignment','center');
set(hcby,'rotation',270,'pos',[3.2 -2.5 0]);


%% SAVE
if ifsave
    fprintf('   saving fig\n');
    save2pdf(25,ofile,'/');
end

end
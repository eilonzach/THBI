function plot_SUITE_of_MODELS( suite_of_models,posterior,ifsave,ofile,lalo )
% plot_SUITE_of_MODELS( suite_of_models,ifsave,ofile )
%
% function to plot final model (and how it compares to the true model)

if nargin<3 || isempty(ifsave)
    ifsave = false; % default is not to save
end
if nargin<4 || isempty(ofile)
    ofile = 'figs/suite_of_models_fig.pdf';
end

sm = suite_of_models;

%% PLOT SUITE OF ACCEPTED MODEL
figure(85); clf; set(gcf,'pos',[120 151 920 947])
ax1 = subplot(1,7,[1:3]); hold on
ax2 = subplot(1,7,[4:6]); hold on
ax4 = subplot(1,7,7); hold on

ichains = unique([sm.chain]);
nchains = length(ichains);

for iii = 1:nchains
ic = find(sm.chain==iii);
if isempty(ic), continue; end
if length(ic)>=200, ic = ic(randperm(length(ic),200)); end;
basecol = colour_get(iii,nchains+1,0,parula);
basecol = 1 - 0.5*(1-basecol); basecol = basecol(:)';

plot(ax1,sm.VS(:,ic),sm.Z,'--','Linewidth',0.05,'color',basecol); hold on
% plot(ax1,PRmodel.vs,sm.Z,'-r','Linewidth',1.5);
set(ax1,'ydir','reverse','fontsize',14,'ytick',[0:25:max(sm.Z)],'color','none');
title(ax1,'Vs','fontsize',20)
xlabel(ax1,'Vs (km/s)','fontsize',16)
ylabel(ax1,'Depth (km)','fontsize',16)
ylim(ax1,[0 max(sm.Z)])

plot(ax2,sm.VP(:,ic),sm.Z,'--','Linewidth',0.05,'color',basecol); hold on
% plot(ax2,PRmodel.vp,sm.Z,'-r','Linewidth',1.5);
set(ax2,'ydir','reverse','fontsize',14,'yticklabel','','ytick',[0:25:max(sm.Z)],'color','none');
title(ax2,'Vp','fontsize',20)
xlabel(ax2,'Vp (km/s)','fontsize',16)
ylim(ax2,[0 max(sm.Z)])
% ylabel(ax2,'Depth (km)','fontsize',16)


Zdh = midpts([0:1:sm.Z(end)]');
nm = hist(posterior.zmoh(ic),Zdh);
ns = hist(posterior.zsed(ic),Zdh);
plot(ax4,nm/sum(nm),Zdh,'k')
fill(ax4,[0,nm/sum(nm),0],[0;Zdh;sm.Z(end)],basecol,'facealpha',0.4)
plot(ax4,ns/sum(ns),Zdh,'k')
fill(ax4,[0,ns/sum(ns),0],[0;Zdh;sm.Z(end)],basecol,'facealpha',0.4)
set(ax4,'ydir','reverse','fontsize',14,'yticklabel','','ytick',[0:25:max(sm.Z)],'color','none');
title(ax4,'Disc.','fontsize',20)
xlabel(ax4,'fraction','fontsize',16)
ylim(ax4,[0 max(sm.Z)])
xmx(iii) = 1.05*max(nm/sum(nm));
end

xlim(ax4,[0,max(xmx)])


%% GET TARGET MODEL for comparison
if nargin > 4 && ~isempty(lalo)
% read 1D profile from Shen and Ritzwoller model
    ncfile = '~/Work/data/models_seismic/SR16_3d_Vs/US.2016.nc';
    [ model ] = read_seismodel_nc( ncfile );

    vs = squeeze(model.Vsv(mindex(model.lon,mod(lalo(2),360)),mindex(model.lat,lalo(1)),:));
    vp = squeeze(model.Vpv(mindex(model.lon,mod(lalo(2),360)),mindex(model.lat,lalo(1)),:));
    Z = model.Z;

else 
    % read in prem
    PRmodel = prem('depths',sm.Z);
    vs = PRmodel.vs;
    vp = PRmodel.vp;
    Z = sm.Z;
end
plot(ax1,vs,Z,'-r','Linewidth',1.5);
plot(ax2,vp,Z,'-r','Linewidth',1.5);

%% SAVE
if ifsave
    fprintf('   saving fig\n')
    save2pdf(85,ofile,'/');
end

end


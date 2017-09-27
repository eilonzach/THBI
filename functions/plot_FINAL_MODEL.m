function plot_FINAL_MODEL( final_model,posterior,ifsave,ofile,ifcomptru,lalo )
% plot_FINAL_MODEL( final_model,posterior,ifsave,ofile )
%
% function to plot final model (and how it compares to the true model)

if nargin<3 || isempty(ifsave)
    ifsave = 0; % default is not to save
end
if nargin<4 || isempty(ofile)
    ofile = [pwd,'/figs/final_model_fig.pdf'];
end
if nargin<5 || isempty(ifcomptru)
    ifcomptru=false;
end

Z = final_model.Z;

%% PLOT FINAL MODEL
figure(95); clf; set(gcf,'pos',[120 308 1100 780])
ax(1) = subplot(1,7,[1,2]); hold on
ax(2) = subplot(1,7,[3,4]); hold on
ax(3) = subplot(1,7,[5,6]); hold on
ax4 = subplot(1,7,7); hold on

variable = {'VS','VP','rho'};
xlims = [[3.2,4.9];[5.9 9.1];[2.6 3.7]];
for iv = 1:length(variable)

fill(ax(iv),[final_model.([variable{iv},'sig1'])(:,1);...
              flipud(final_model.([variable{iv},'sig1'])(:,2))],...
              [Z;flipud(Z)],'-','Linewidth',1.5,'Facecolor',[0.7 0.7 0.7],'Edgecolor',[0.6 0.6 0.6]);
plot(ax(iv),final_model.([variable{iv},'av']),Z,'-r','Linewidth',2);
plot(ax(iv),final_model.([variable{iv},'sig2']),Z,'-','color',[0.4 0.4 0.4],'Linewidth',1);
set(ax(iv),'ydir','reverse','fontsize',14,'ytick',[0:25:max(Z)],'ylim',[0 max(Z)],'xlim',xlims(iv,:),...
    'color','none','box','on');
title(ax(iv),variable{iv},'fontsize',20)
xlabel(ax(iv),[variable{iv},' (km/s)'],'fontsize',16)

if iv==1
    ylabel(ax(iv),'Depth (km)','fontsize',16) 
    % moh
    Zmoh(1) = final_model.Zd(2).mu;
    Zmoh(2) = final_model.Zd(2).std;
    plot(ax(iv),xlims(1,:),Zmoh(1)*[1 1],'--k','linewidth',2)
    text(ax(iv),xlims(1,1)+0.1,Zmoh(1)+8,sprintf('$%.1f \\pm %.1f$',Zmoh),...
        'fontsize',17,'interpreter','latex')    
else
    set(ax(iv),'yticklabel','')
end

end

%% histogram of discontinuities
Zdh = midpts([0:0.5:max(Z)]); %Zdh = Zdh(1:3:end);
%moho
nm = hist(posterior.zmoh,Zdh);
plot(ax4,nm/sum(nm),Zdh,'k')
fill(ax4,nm/sum(nm),Zdh,[0.5 0.5 0.5])
% seds
if ~all(posterior.zsed==0)
ns = hist(posterior.zsed,Zdh);
plot(ax4,ns/sum(ns),Zdh,'k')
fill(ax4,[0,ns/sum(ns)],[0,Zdh],[0.5 0.5 0.5])
end
set(ax4,'ydir','reverse','fontsize',14,'yticklabel','','ytick',[0:25:max(Z)],'ylim',[0 max(Z)],'color','none');
title(ax4,'Disc.','fontsize',20)
xlabel(ax4,'fraction','fontsize',16)
xlim(ax4,[0 max(nm/sum(nm))])

%% GET TARGET MODEL for comparison
if ifcomptru
if nargin > 5 && ~isempty(lalo)
% read 1D profile from Shen and Ritzwoller model
    ncfile = '~/Work/data/models_seismic/SR16_3d_Vs/US.2016.nc';
    [ model ] = read_seismodel_nc( ncfile );

    vs = squeeze(model.Vsv(mindex(model.lon,mod(lalo(2),360)),mindex(model.lat,lalo(1)),:));
    vp = squeeze(model.Vpv(mindex(model.lon,mod(lalo(2),360)),mindex(model.lat,lalo(1)),:));
    Z = model.Z;

else 
    % read in prem
    PRmodel = prem('depths',Z);
    vs = PRmodel.vs;
    vp = PRmodel.vp;
end
plot(ax(1),vs,Z,'-b','Linewidth',1.5);
plot(ax(2),vp,Z,'-b','Linewidth',1.5);
end

%% SAVE
if ifsave
    save2pdf(95,ofile,'/');
end

end


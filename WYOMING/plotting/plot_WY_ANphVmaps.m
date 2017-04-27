%% Simple script to plot the ambient noise phase velocity maps 
% from Shen and Ritzwoller 2016

latlims = [40 48.8];
lonlims = [-113 -103];

addpath('~/Documents/MATLAB/BayesianJointInv/WYOMING/matguts');

datadir = '~/Work/data/models_seismic/US_RAYLEIGH_PHASE_ANT_VEL_SHEN/'; % need final slash

periods = [8:2:32,36,40]'; freqs = 1./periods;
Nf = length(freqs);

% specify frequencies to actually plot, to match RD16 supplemental fig. 4
periods2plot = [8:4:36];

phV_mean = nan(Nf,1);
figure(1), clf; set(gcf,'pos',[700 125 740 900])
for ii = 1:Nf
   
    [ lat,lon,phV, latgrid, longrid,phVgrid ] = load_ANphV_data( periods(ii),datadir );

    % average phase velocity at this freq. in well resolved region
    phV_mean(ii) = nanmean(nanmean(phVgrid));
    
    % phase velocity deviation from mean at that freq., in percent
    dphVgrid = 100*(phVgrid-phV_mean(ii))./phV_mean(ii); 
    
    if any(periods2plot==periods(ii))

        figure(1), 
        ax = subplot(3,3,9-find(periods2plot==periods(ii)));
        contourf(ax,longrid,latgrid,dphVgrid,20,'color','none')
        shading flat
        caxis([-6 4])
        colormap(flipud(jet))
        set(ax,'ylim',latlims,'xlim',lonlims,'fontsize',12)
        add_state_boundaries(ax,latlims,lonlims)
    end    
    
    
end

% dispersion curves
subplot(3,3,9), cla,hold on
% average
plot(periods,phV_mean,'k-o'); axis([0 40 2.9 4.1])

% yellowstone point
[ iphV ] = disp_curve_AN_latlon( periods,45,-111 );
plot(periods,iphV,'r-o')
% craton point
[ iphV ] = disp_curve_AN_latlon( periods,44,-106 );
plot(periods,iphV,'b-o')% 

% colorbar
colorbar('location','southoutside')
caxis([-6 4])
colormap(flipud(jet))


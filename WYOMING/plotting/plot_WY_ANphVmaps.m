%% Simple script to plot the phase velocity maps 
% from Riddhi Dave and Aibing Li (Geology 2016)

latlims = [40 50];
lonlims = [-113 -103];

addpath('~/Documents/MATLAB/BayesianJointInv/WYOMING/matguts');

datadir = '~/Work/data/models_seismic/US_RAYLEIGH_PHASE_ANT_VEL_SHEN/'; % need final slash

periods = [8:2:40]'; freqs = 1./periods;
Nf = length(freqs);

% specify frequencies to actually plot, to match RD16 supplemental fig. 4
periods2plot = [8 16 40];

phV_mean = nan(Nf,1);
figure(1), clf; set(gcf,'pos',[700 125 740 900])
for ii = 1:Nf
    datafile = [num2str(periods(ii)),'_ANT.vel'];
    dat = load([datadir,datafile]);
    lon = dat(:,1); lat = dat(:,2); phV = dat(:,3);
    scatter(lon,lat,40,phV)
    return
    [ lat,lon,phV, latgrid, longrid,phVgrid ] = load_phV_data( freqs(ii),datadir );

    % fractional error, in percent
    frac_err_grid = 100*errgrid./phVgrid; 
    
    % average phase velocity at this freq. in well resolved region
    phV_mean(ii) = mean(mean(phVgrid(frac_err_grid<1.5)));
    
    % phase velocity deviation from mean at that freq., in percent
    dphVgrid = 100*(phVgrid-phV_mean(ii))./phV_mean(ii); 
    
    if any(periods2plot==freqs(ii))
        % mask hi-error region
        dphVgrid(frac_err_grid20>1.5) = nan;

        figure(1), 
        ax = subplot(3,3,9-find(periods2plot==freqs(ii)));
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
plot(periods,phV_mean,'k-o'); axis([0 170 3.4 4.6])

% yellowstone point
[ iphV ] = disp_curve_latlon( freqs,45,-111 );
plot(periods,iphV,'r-o')
% craton point
[ iphV ] = disp_curve_latlon( freqs,44,-106 );
plot(periods,iphV,'b-o')% 

% colorbar
colorbar('location','southoutside')
caxis([-6 4])
colormap(flipud(jet))


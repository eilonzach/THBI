%% Simple script to plot the phase velocity maps 
% from Riddhi Dave and Aibing Li (Geology 2016)

latlims = [40 50];
lonlims = [-113 -103];

addpath('~/Documents/MATLAB/BayesianJointInv/WYOMING/matguts');

datadir = '~/Dropbox/Dave_Li_phV/2D_phase_velocities/'; % need final slash
errdir = '~/Dropbox/Dave_Li_phV/Errors/'; % need final slash

[ freqs, periods ] = get_freqs(datadir);
Nf = length(freqs);

% specify frequencies to actually plot, to match RD16 supplemental fig. 4
freqs2plot = [.006 0.008 0.01 0.015 0.02 0.029 0.037 0.05 ];

% get 20s error map for the masking - as RD16 did (sneakily)
[ ~,~,~,~,~,phVgrid ] = load_phV_data( 1./20,datadir );
[ ~,~,~,~,~,errgrid ] = load_phV_error( 1./20,errdir );
frac_err_grid20 = 100*errgrid./phVgrid; 

phV_mean = nan(Nf,1);
figure(1), clf; set(gcf,'pos',[700 125 740 900])
for ii = 1:Nf
    [ lat,lon,phV, latgrid, longrid,phVgrid ] = load_phV_data( freqs(ii),datadir );
    [ ~,~,err, ~, ~,errgrid ] = load_phV_error( freqs(ii),errdir );

    % fractional error, in percent
    frac_err_grid = 100*errgrid./phVgrid; 
    
    % average phase velocity at this freq. in well resolved region
    phV_mean(ii) = mean(mean(phVgrid(frac_err_grid<1.5)));
    
    % phase velocity deviation from mean at that freq., in percent
    dphVgrid = 100*(phVgrid-phV_mean(ii))./phV_mean(ii); 
    
    if any(freqs2plot==freqs(ii))
        % mask hi-error region
        dphVgrid(frac_err_grid20>1.5) = nan;

        figure(1), 
        ax = subplot(3,3,9-find(freqs2plot==freqs(ii)));
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


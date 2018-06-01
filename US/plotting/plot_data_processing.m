close all 
clear 
%% setup
datadir = './'; % need final slash
station = 'RSSD';
network = 'IU';
gc_av = 71; % 
seaz_av = 148; % 
ifsave = true;


%% colours
cl1 = [0.8 0.8 1];
cl2 = [1 0.8 0.8];
cols = [cl1;cl2];

%% process levels
proclevs = [2,3,4,6];
ylabs = {'Raw',...
         'Stage 1',...
         'Stage 2',...
         'Final'};

%% load data
load(sprintf('AVARS/plotting_dat_%s_%s_%.0f_%.0f',station,network,gc_av,seaz_av));

phases = plot_dat{1,3};
components = plot_dat{2,3};

%% fig setup
figure(5); clf, set(gcf,'pos',[150 13 1200 792])
for ip = 1:length(phases)
    for iproc = 1:4
        xl = 0.05 + 0.47*(ip-1);
        yb = 0.08 + 0.22*(4-iproc);
        ax(ip,iproc) = axes(gcf,'pos',[xl yb 0.43 0.19]); hold on %#ok<SAGROW>
    end
end

%% get pol to make pos. spike
pold = zeros(length(phases),1);
pols = zeros(length(phases),1);
for ip = 1:length(phases)
    pold(ip) = sign(maxab(maxab(squeeze(sum(plot_dat{6,ip}.data,2)))'));
    pols(ip) = sign(maxab(maxab(plot_dat{8,ip}.data)'));
end

%% plot the background waveforms
for ip = 1:length(phases)
    for iproc = 1:4
        proclev = proclevs(iproc);
        pol = sign(maxab(squeeze(maxab(plot_dat{proclev,ip}.data))'));
        for ic = 1:2
            % plot the data!
            if iproc~=4
            plot(ax(ip,iproc),...
                plot_dat{proclev,ip}.tt(:,ic),...
                pold(ip)*plot_dat{proclev,ip}.data(:,:,ic),'color',cols(ic,:));
            else
            % use aligned time axis if last one
            plot(ax(ip,iproc),...
                plot_dat{8,ip}.tt,...
                pold(ip)*plot_dat{proclev,ip}.data(:,:,ic),'color',cols(ic,:));
            end
        end
    end
end

%% plot the stack waveforms
for ip = 1:length(phases)
for ic = 1:2
plot(ax(ip,4),...
     plot_dat{8,ip}.tt,...
     pols(ip)*plot_dat{8,ip}.data(:,ic),'color',cols(ic,:).^4, 'linewidth',2);
end
end

%% limits and things
for ip = 1:length(phases)
    for iproc = 1:4
        % axis limits, fontsize
        set(ax(ip,iproc),'fontsize',15,'box','on','linewidth',1.3,...
            'xlim',50*[-1 1],'xtick',[-100:10:100],...
            'ylim',1.1*[-1 1],'ytick',[-1:1:1],...
            'layer','top')
        
        % yticks, ylabel
        if ip == 2
            set(ax(ip,iproc),'yticklabel',[])
        else
            ylabel(ax(ip,iproc),ylabs{iproc},...
                'interpreter','latex','fontsize',22)            
        end
        
        % xticks,xlabel
        if iproc ~= 4
            set(ax(ip,iproc),'xticklabel',[])
        else
            xlabel(ax(ip,iproc),sprintf('Time from %s arrival',phases{ip}(1)),...
                'interpreter','latex','fontsize',22)
        end
 
        % Title
        if iproc == 1   
            title(ax(ip,iproc),sprintf('\\textbf{%s data}',phases{ip}),...
                'interpreter','latex','fontsize',24)            
        end
    end
end

%% N traces
for ip = 1:length(phases)
    for iproc = 1:4
        proclev = proclevs(iproc);
            % report number of traces
            figN_add([num2str(size(plot_dat{proclev,ip}.data(:,:,1),2)) ' traces'],...
                     ax(ip,iproc),0.03,0.84,12,'none','bold');
    end
end
   

%% Legend
axl = axis(ax(1,4));
xlr = axl(1) + [0.04,0.12]*diff(axl(1:2));
ybt = axl(3) + [0.13,0.21]*diff(axl(3:4));
for ic = 1:2
    plot(ax(1,4),xlr,ybt(ic)*[1 1],'color',cols(ic,:).^4,'linewidth',2)
    text(ax(1,4),xlr(end)+0.1*diff(xlr),ybt(ic),components{ic},'fontweight','bold','fontsize',13)
end

%% Station, gcarc etc info
if length(phases)==2
    text(ax(2,4),xlr(1),mean(ybt),...
     sprintf('%s \\,  $\\Delta \\sim$%.0f$^{\\circ}$  \\,  $\\phi \\sim$%.0f$^{\\circ}$',station,gc_av,seaz_av),...
     'interpreter','latex','fontsize',13)
else
    text(ax(1,4),mean(xlr)+1.5*diff(xlr),mean(ybt),...
     sprintf('%s \\,  $\\Delta \\sim$%.0f$^{\\circ}$  \\,  $\\phi \\sim$%.0f$^{\\circ}$',station,gc_av,seaz_av),...
     'interpreter','latex','fontsize',13)
end
 
%% Save
if ifsave
    save2jpg(5,'data_processing')
end






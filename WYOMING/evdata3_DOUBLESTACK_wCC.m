clear all
close all
station = 'RSSD';
network = 'IU';
gc_lims = [65,75];

RFphase = {'Ps','Sp'};

xcor_filtfs = [[0.1;2],[0.05;1]];
% xcor_filtfs = [[0.02;100],[0.02;100]];
xcor_win = [-20 10];
npoles = 2;

filtfs = [0.02;100]; % filter frequencies pre PCA, in Hz

datwind = [-50 50]; % in sec around main arrival
tapertime = 5; % s at beginnning and end of window to taper over
acormin = 0.8;

SNRmin = 5;
maxdepth = 50;  % maximum EQ depth
min4clust = 20; % minimum # of EQ to even try the stack
min4stack = 5;  % minumum # of traces in order to keep the stack

ifsave = false;
verbose = true;

load(sprintf('dat_%s_%s_%.0fto%.0f',station,network,gc_lims(1),gc_lims(2)));
addpath('matguts')

%% cluster analysis on events
Z = linkage([eqar.gcarcs.*sind(eqar.seazs),eqar.gcarcs.*cosd(eqar.seazs)]);
T = cluster(Z,'cutoff',5,'criterion','distance');
nT = histcounts(T);
iclusts = find(nT>min4clust);

% plot EQ and clusters
addpath('~/Documents/MATLAB/lib/m_map')
figure(78), clf, hold on, set(gcf,'color','w','pos',[1360 685 560 420]) 
m_proj('ortho','lat',eqar.slat,'long',eqar.slon);
m_grid('xtick',[],'ytick',[],'linestyle','-','backcolor',[.9 .99 1]);
m_coast('patch',[255,225,204]/255,'edgecolor','none');
for icl = 1:length(iclusts)
    ind = (T==iclusts(icl));
    m_line([eqar.slon,mean(eqar.elons(ind))],[eqar.slat,mean(eqar.elats(ind))],'color','r','linewi',2);
    m_plot(eqar.elons(ind),eqar.elats(ind),'o','markersize',10,'markerfacecolor',mean(colour_get(T(ind),length(nT),1,parula)),'markeredgecolor','none')
end
m_scatter(eqar.elons,eqar.elats,70,colour_get(T,length(nT),1,parula));
m_line(eqar.slon,eqar.slat,'marker','^','color','k','markerfacecolor','k','markersize',10);
m_range_ring(eqar.slon,eqar.slat,[1100:1100:10000],'color','k','linewi',1);
xlabel(['10 deg range rings from ',eqar.sta],'fontsize',20);

save2pdf(78,['EQ_source_',eqar.sta],'figs')

%% loop through clusters
for icl = 1:length(iclusts)

iclust = iclusts(icl);

for ip = 1:length(RFphase)
    fprintf('\nCluster %.0f, %s \n',icl,RFphase{ip})
    switch RFphase{ip}(1)
        case 'P', chp = 1; chd = 2; chidx = 1;% main phase is Z, align on Z
        case 'S', chp = 2; chd = 1; chidx = 2;% main phase is R, align on T
    end
    
%% =================  ALIGN ON CHIDX COMPONENT USING CC  ==================
    evinds = 1:eqar.norids;
    % QC on depth, only include this cluster
    edeps = [eqar.edeps];
	evinds(edeps(evinds)>=maxdepth) = [];
	evinds(T(evinds)~=iclust) = [];
    plot_rec_sec_align(eqar,RFphase{ip}(1),eqar.components{chp},filtfs,evinds)
    
    %% Make matrix of data to cross-correlate        
    % grab the sample rate
    samprate = unique(round(1./diff(eqar.tt(:,:,ip))));
    % filter before cross correlation
    xcordat = eqar.dataZRT(:,:,chidx,ip);
    xcordat = filt_quick(xcordat,xcor_filtfs(1,ip),xcor_filtfs(2,ip),1./samprate,npoles);
    % taper extreme ends
%     xcordat = flat_hanning_win(xcortt,xcordat,xcortt(1)+5,xcortt(end)-5,5);
    % detrend
    xcordat = detrend(xcordat);
    % timing
    xcortt = eqar.tt(:,1,ip);
    xcinds = xcortt>=xcor_win(1) & xcortt<xcor_win(2); % indices of times within window
    % taper
    xcortt = eqar.tt(:,1,ip);
    xcordat = flat_hanning_win(xcortt,xcordat,xcortt(1)+5,xcortt(end)-5,5);
    
    % SNR after the filter
    nswin = find(xcortt>(max([datwind(1)+tapertime,-55])) & xcortt<=-5); % noise window 
    dwin  = find(xcortt>=-5 & xcortt<=10); % generous first arrival window 
    SNR_est = max(abs(xcordat(dwin,:)))./2./rms(xcordat(nswin,:));

    % window to xcorr time segment
    xcordat = detrend(xcordat(xcinds,:),'constant');
    
    % normalise
    for ie = 1:eqar.norids,
        dd = xcordat(:,ie);
        ddmax = max(abs(dd));
        xcordat(:,ie) = dd./ddmax; % flip and scale 
%         xcordat(:,ie) = dd.*ddsgn; % flip based on sign
%         xcordat(:,ie) = dd./maxab(dd); % flip and scale 
    end
    % window before cross correlation to xcor_window around arrival
    xcordat = flat_hanning_win(xcortt(xcinds),xcordat,xcor_win(1),xcor_win(2),1);

    %% polarity
    indtemplate = evinds(mindex(-SNR_est(evinds)));
    ddpol = polarity_est(xcordat(:,indtemplate),xcordat(:,evinds),5*samprate);
    ddpol(evinds) = ddpol; if length(ddpol)<eqar.norids, ddpol(eqar.norids)=0; end
    % flip depending on polarity of maximum
    for ie = 1:eqar.norids,
        xcordat(:,ie) = xcordat(:,ie).*ddpol(ie);
    end
    evinds(ddpol(evinds)==0) = [];

    %% CROSS CORRELATE - iteratively throw out low acor
    acor = -1;
    evinds(SNR_est(evinds)<2) = [];
    
    figure(2),clf;plot(xcortt(xcinds),xcordat(:,evinds))

    while any(acor<0)
        fprintf('xcorr-ing %.0f traces\n',length(evinds))
        [dcor, dcstd, dvcstd, acor]=xcortimes(xcordat(:,evinds),1./samprate, -xcor_win(1), 5,verbose);
        fprintf('killing %.0f traces\n',sum(acor<0))
        evinds(acor<0) = [];
    end
    while any(acor<acormin)
        fprintf('xcorr-ing %.0f traces\n',length(evinds))
        
        [dcor, dcstd, dvcstd, acor]=xcortimes(xcordat(:,evinds),1./samprate, -xcor_win(1), 5,verbose);
        fprintf('killing %.0f traces\n',sum(acor<acormin))
        evinds(acor<acormin) = [];
    end
    
    plot_rec_sec_align(eqar,RFphase{ip}(1),eqar.components{chp},filtfs,evinds,-dcor,ddpol)
%     plot_rec_sec_align(eqar,RFphase{ip}(1),eqar.components{chp},xcor_filtfs,evinds,-dcor,ddpol)
    
    % subset to only good events
    ddpol = ddpol(evinds);
    
    dcor_all = nan(eqar.norids,1);
    dcor_all(evinds) = dcor;

% =============================  END CC  ================================
    
    
%% PROCESS data traces: shift by dcor, filter, window, normalise by max on chid
    
    dataZRT_proc = nan([size(eqar.dataZRT(:,:,1,ip)),3]);
    SNR_all = zeros(length(evinds),3);
    for ie = 1:length(evinds)
        evind = evinds(ie);
        % get data for this evt
        dd = squeeze(eqar.dataZRT(:,evind,:,ip));
        % flip sign if needed 
%         dd = dd.*ddpol(ie); % flip sign if needed
        % shift according to cross correlation 
        dd_s = interp1(eqar.tt(:,evind,ip), dd ,eqar.tt(:,evind,ip)+dcor(ie)); % interp to new time axis with the shift
        % set any extrapolated vals (made nan) to zero
        dd_s(isnan(dd_s)) = 0;
        % Filter:
        dd_sf = filt_quick(dd_s,filtfs(1),filtfs(2),1./samprate);
        % taper 
        dd_sft = flat_hanning_win(eqar.tt(:,evind,ip),dd_sf,datwind(1),datwind(2),tapertime);
        % nominal arrival window

        normval = max(max(abs(dd_sft)));
        dd_sftn   = dd_sft./normval./ddpol(ie); % scale and flip sign

        % plot
%         figure(1),subplot(511),plot(dd),subplot(512),plot(dd_s),subplot(513),plot(dd_sf),subplot(514),plot(dd_sft),subplot(515),plot(dd_sftn)

        % insert into _shifted, _normalised structure
        dataZRT_proc(:,evind,:,ip) = dd_sftn;
        
        %% SNR
        nswin = find(eqar.tt(:,evind,ip)<=-5 & eqar.tt(:,evind,ip)>(max([datwind(1)+tapertime,-55]))); % noise window 
        dwin  = find(eqar.tt(:,evind,ip)>=-2 & eqar.tt(:,evind,ip)<=10); % generous first arrival window 
        for ic = 1:3
            SNR_all(ie,ic) = max(abs(detrend(dd_sftn(dwin,ic))))./2./rms(detrend(dd_sftn(nswin,ic)));
        end
    end
    
    % kill bad SNR
    evinds(SNR_all(:,chp)<SNRmin | SNR_all(:,chd)<SNRmin/10) = [];
    

    %% FIRST-STACK trace
    for ic = 1:3
        dataZRT_stk(:,ic,ip) = sum(dataZRT_proc(:,evinds,ic,ip),2)/length(evinds); %#ok<SAGROW>
    end
    
    %% Xcorr again!
    for ie = 1:length(evinds)
        evind = evinds(ie);        
        dd = squeeze(eqar.dataZRT(:,evind,:,ip));
        % fix
        dd_s = interp1(eqar.tt(:,evind,ip), dd ,eqar.tt(:,evind,ip)+dcor(ie)); % interp to new time axis with the shift
        dd_s(isnan(dd_s)) = 0;
        % filt
        dd_sf = filt_quick(dd_s,filtfs(1),filtfs(2),1./samprate);
        % taper
        dd_sft = flat_hanning_win(eqar.tt(:,evind,ip),dd_sf,datwind(1),datwind(2),tapertime);
        normval = max(max(abs(dd_sft)));
        pol = polarity_est(dataZRT_stk(:,chp,ip),...
                           flat_hanning_win(eqar.tt(:,evind,ip),dd_sft(:,chp),-20,20,tapertime),...
                            5*samprate);
        dd_sftnp = dd_sft./normval./pol;
        pause;
    end

    

    %% pick phase onset
    figure(20), clf, set(gcf,'pos',[15 25 1426 420]), hold on
    htr = plot(eqar.tt(:,1,ip),dataZRT_stk(:,:,ip),'linewidth',1.5);
    title('PICK PHASE ONSET','fontsize',22)
    set(gca,'xlim',[-50 50],'ylim',[-1. 1.],'fontsize',18)
    legend(eqar.components)
    tshft = ginput(1);
    tt_stk(1:size(eqar.tt,1),ip) = eqar.tt(:,1,ip) - tshft(1);
    tt_stk = round_level(tt_stk,1/samprate);
    for ic = 1:length(htr)
        htr(ic).XData = tt_stk(:,ip);
        htr(ic).YData = dataZRT_stk(:,ic,ip);
    end
    pause(0.1)

    %% other averaged parms
    gcarc_av(ip) = mean(eqar.gcarcs(evinds));
    seaz_av(ip) = mean(eqar.seazs(evinds));
    rayp_av(ip) = mean(eqar.rayps(evinds,ip));
    norids(ip,1) = length(evinds);
    evinds_save{ip,1} = evinds;
    
    if length(evinds) < min4stack, continue; end
    
    %% nice plot
    figure(21), clf, set(gcf,'pos',[92 404 802 642])
    for ic = 1:length(eqar.components)
        subplot(3,1,ic), hold on
        plot(tt_stk(:,ip),dataZRT_proc(:,evinds,ic,ip),'color',[.8 .89 .9],'linewidth',0.5)
        plot(tt_stk(:,ip),dataZRT_stk(:,ic,ip),'k','linewidth',3)
        set(gca,'xlim',[-50 50],'ylim',[-1.1 1.1],'fontsize',16)
        ylabel([eqar.components{ic},'-comp'],'fontsize',20)
    end
    subplot(311),title([RFphase{ip},' stack'],'fontsize',22)
    subplot(313),xlabel('Time (s)','fontsize',20)
    if ifsave
        save2jpg(21,sprintf('stackdata_%s_%s_%.0f_%.0f_%s',eqar.sta,eqar.nwk,gcarc_av(ip),seaz_av(ip),RFphase{ip}),'figs')
    end
    pause
end % loop on phases

%% kill ones with too few evts in each stack
keepstk = norids >= min4stack;
evinds_save = evinds_save(keepstk);
evinds_save = evinds_save(keepstk);
seaz_av = seaz_av(keepstk);
gcarc_av = gcarc_av(keepstk);
rayp_av = rayp_av(keepstk);
dataZRT_stk = dataZRT_stk(:,:,keepstk);
tt_stk = tt_stk(:,keepstk);
phases = {eqar.phases(keepstk)};

%% save
avar = struct('sta',eqar.sta,'nwk',eqar.nwk,'phases',phases,'components',{eqar.components},...
              'slat',eqar.slat,'slon',eqar.slon,'selev',eqar.selev,...
              'norids',norids,'evinds',{evinds_save},...
              'elats',eqar.elats,'elons',eqar.elons,'edeps',eqar.edeps,...
              'emags',eqar.emags,'evtimes',{eqar.evtimes},...
              'seaz',seaz_av,'seazs',eqar.seazs,...
              'gcarc',gcarc_av,'gcarcs',eqar.gcarcs,...
              'rayp',rayp_av,'rayps',eqar.rayps,...
              'SNRmin',SNRmin,'acormin',acormin,...
              'xcor_filt',xcor_filtfs,'xcor_win',xcor_win,...
              'filtfs_PCA',filtfs,...
              'dataZRT',dataZRT_stk,'tt',tt_stk);
          
if ifsave
    fprintf('SAVING\n')
    save(sprintf('avar_dat_%s_%s_%.0f_%.0f',station,network,mean(gcarc_av),mean(seaz_av)),'avar');
else
    pause 
end

end % loop on clusters
    

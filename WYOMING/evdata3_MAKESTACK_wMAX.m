% Align on the maximum of the arrival waveform for each event. 

clear all
close all
station = 'LKWY';
network = 'US';
gc_lims = [65,75];

RFphase = {'Ps','Sp'};

xcor_filtfs = [[0.1;2],[0.05;1]];
xcor_win = [-20 8];
npoles = 4;

filtfs = [0.02;100]; % filter frequencies pre PCA, in Hz

datwind = [-50 50]; % in sec around main arrival
tapertime = 5; % s at beginnning and end of window to taper over
acormin = 0.8;

SNRmin = 5;
maxdepth = 100;

min4clust = 20;

ifsave = true;
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
figure(78), clf, hold on, set(gcf,'color','w') 
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

for icl = 1:length(iclusts)

iclust = iclusts(icl);

for ip = 1:length(RFphase)
    switch RFphase{ip}(1)
        case 'P', chp = 1; chd = 2; chidx = 1;% parent phase is Z, align on Z
        case 'S', chp = 2; chd = 1; chidx = 2;% main phase is R, align on T
    end
    
%% =================  ALIGN ON CHIDX COMPONENT USING CC  ==================
    evinds = 1:eqar.norids;
    % QC on depth, only include this cluster
    edeps = [eqar.edeps];
	evinds(edeps(evinds)>=maxdepth) = [];
	evinds(T(evinds)~=iclust) = [];

    
    % Make matrix of data to find arrival times        
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
    figure(1),clf;plot(xcortt,xcordat)
    
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

    % polarity
    indtemplate = evinds(mindex(-SNR_est(evinds)));
    ddpol = polarity_est(xcordat(:,indtemplate),xcordat,5*samprate);
    % flip depending on polarity of maximum
    for ie = 1:eqar.norids,
        xcordat(:,ie) = xcordat(:,ie).*ddpol(ie);
    end
    evinds(ddpol(evinds)==0) = [];
    
    % CROSS CORRELATE - iteratively throw out low acor
    evinds(SNR_est(evinds)<2) = [];
    
    wintt = xcortt(xcinds);
    dcor = zeros(length(evinds),1);
    for ie = 1:length(evinds)
        dcor(ie) = wintt(mindex(xcordat(:,evinds(ie)),1));
    end
	dcor = dcor-mean(dcor);
    % kill too-large dcors
    evinds(abs(dcor)>3) = [];
    dcor(abs(dcor)>3) = [];
    
    plot_rec_sec_align(eqar,'P','Z',filtfs,evinds,-dcor,ddpol)
    plot_rec_sec_align(eqar,'P','Z',xcor_filtfs,evinds,-dcor,ddpol)
    
    
    % subset to only good events
    ddpol = ddpol(evinds);
    
    dcor_all = nan(eqar.norids,1);
    dcor_all(evinds) = dcor;
%     acor_all = nan(eqar.norids,1);
%     acor_all(evinds) = acor;
%     
% =============================  END CC  ================================
    
    
%% PROCESS data traces: shift by dcor, filter, window, normalise by max on chid
    
    dataZRT_proc = nan([size(eqar.dataZRT(:,:,1,ip)),3]);
    SNR_all = zeros(length(evinds),3);
    figure(99);
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
%         dwin = find(eqar.tt(:,evind,ip)>=-10 & eqar.tt(:,evind,ip)<=5);
        % normalise by the maximum value on the parent component
%         Imx = maxab(dd_sft(dwin,chid));
%         Xmx = maxab(dd_sft(dwin,chidx));
%         normval = max(abs(dd_sft(:,chid)));
%         dd_sftn(:,1:2) = dd_sft(:,1:2)./normval;  % scale and flip sign based on max Z or R
%         dd_sftn(:,3)   = dd_sft(:,3)  ./normval; % scale and flip sign  based on max Z or T
        dd_sftn   = dd_sft./max(abs(dd_sft(:,chp)))./ddpol(ie); % scale and flip sign
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
        
%         [Imx,
        if SNR_all(ie,chp)>=SNRmin
            for icc =1:3,subplot(3,1,icc),hold on,plot(dd_sftn(:,icc));end
        end

    end
    
    evinds(SNR_all(:,chp)<SNRmin) = [];
    evinds(any(abs(dataZRT_proc(:,evinds,chd))>1)) = [];
    
    % plot shifted, normalised data
    if verbose
    figure(6);plot(eqar.tt(:,evinds,ip),dataZRT_proc(:,evinds,chp,ip))
    figure(7);plot(eqar.tt(:,evinds,ip),filt_quick(dataZRT_proc(:,evinds,chp,ip),xcor_filtfs(1,ip),xcor_filtfs(2,ip),1./samprate,npoles))
    figure(8);plot(eqar.tt(:,evinds,ip),filt_quick(dataZRT_proc(:,evinds,chd,ip),xcor_filtfs(1,ip),xcor_filtfs(2,ip),1./samprate,npoles))
    end
    %% STACK trace
    for ic = 1:3
        dataZRT_stk(:,ic,ip) = sum(dataZRT_proc(:,evinds,ic,ip),2)/length(evinds); %#ok<SAGROW>
    end

    if verbose
        figure(9);
        for ic = 1:3
%             subplot(3,1,ic);plot(eqar.tt(:,1,ip),dataZRT_stk(:,ic,ip))
            subplot(3,1,ic);plot(eqar.tt(:,1,ip),filt_quick(dataZRT_stk(:,ic,ip),xcor_filtfs(1,ip),xcor_filtfs(2,ip),1./samprate,npoles))
        end
    end
    
%     %% EOF
%     cvar = zeros(length(evinds),length(evinds),3);
%     PCAmat = zeros(size(dataZRT_proc(:,evinds,:,ip)));
%     pct_var = zeros(length(evinds),3);
%     % break down to PCA
%     for ic = 1:3
%         [cvar(:,:,ic), PCAmat(:,:,ic), ~, ~, pct_var(:,ic)] = pca(dataZRT_proc(:,evinds,ic,ip));
%     end
%     
%     cutoff_mat = diag(cumsum(mean(pct_var,2))<90); % cutoff at val by which have explained 95% of variance
%     cutoff_mat = diag([1,1,1,1,1,zeros(1,length(pct_var)-5)]); % cutoff at val by which have explained 95% of variance
%     N_cutoff = sum(diag(cutoff_mat));
%     
%     % reconstitute
%     dataZRT_procr = zeros(size(dataZRT_proc));
%     for ic = 1:3
%         dataZRT_procr(:,evinds,ic,ip) = PCAmat(:,:,ic)*(cutoff_mat*cvar(:,:,ic))';
%     end    
%     
%     figure(8); clf
%     for ic = 1:3
%         subplot(3,1,ic),plot(eqar.tt(:,1,ip),filt_quick(dataZRT_procr(:,evinds,ic,ip),filtfs(1),filtfs(2),1./samprate,npoles))
%     end
% 
%     %% Stack
%     for ic = 1:3
%         dataZRT_stk(:,ic,ip) = sum(dataZRT_procr(:,evinds,ic,ip),2)/N_cutoff; %#ok<SAGROW>
%     end;
%     figure(8), 
%     for ic = 1:3, 
%         subplot(3,1,ic),hold on
% %         plot(eqar.tt(:,1,ip),filt_quick(dataZRT_stk(:,ic,ip),filtfs(1),filtfs(2),1./samprate,npoles),'k','linewidth',2)
%         plot(eqar.tt(:,1,ip),dataZRT_stk(:,ic,ip),'k','linewidth',2)
%         xlim(datwind)
%     end 
%     
    %% pick phase onset
    figure(20), clf, set(gcf,'pos',[371 678 1426 420])
    plot(eqar.tt(:,1,ip),dataZRT_stk(:,:,ip))
    title('PICK PHASE ONSET','fontsize',22)
    tshft = ginput(1);
    tt_stk(1:size(eqar.tt,1),ip) = eqar.tt(:,1,ip) - tshft(1);
    plot(tt_stk(:,ip),dataZRT_stk(:,:,ip))
    set(gca,'xlim',[-50 50],'ylim',[-1.1 1.1],'fontsize',18)
    legend(eqar.components)

    %% other averaged parms
    gcarc_av(ip) = mean(eqar.gcarcs(evinds));
    seaz_av(ip) = mean(eqar.seazs(evinds));
    rayp_av(ip) = mean(eqar.rayps(evinds,ip));
    norids(ip,1) = length(evinds);
    evinds_save{ip,1} = evinds;
    
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
    save2jpg(21,sprintf('stackdata_%s_%s_%.0f_%.0f_%s',eqar.sta,eqar.nwk,gcarc_av(ip),seaz_av(ip),RFphase{ip}),'figs')
end

%% save
avar = struct('sta',eqar.sta,'nwk',eqar.nwk,'phases',{eqar.phases},'components',{eqar.components},...
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
end

end
    

clear all
close all
datadir = '~/Documents/MATLAB/THBI_paper/DATA/'; % need final slash
station = 'RSSD';
network = 'IU';
gc_lims = [65,75];
% gc_lims = [45,55];

RFphase = {'Ps','Sp'};
% RFphase = {'Ps'};


%% processing parms
npoles = 2;
tapertime = 5; % s at beginnning and end of window to taper over

dpol_win = [-1 5];
datwind = [-50 50]; % in sec around main arrival

filtfs = [0.02;100]; % filter frequencies pre PCA, in Hz

SNRmin = 3;

%% xcor parms
xcor_filtfs = [[0.02;2],[0.05;1]];
% xcor_filtfs = [[0.02;100],[0.02;100]];
xcor_win = [-20 5];
acormin = 0.175;
maxlag = 5;

maxdepth = 50;  % maximum EQ depth
min4clust = 15; % minimum # of EQ to even try the stack
clust_cutoff = 5; % minimum Delta of clustered EQ locations (deg)
min4stack = 5;  % minumum # of traces in order to keep the stack

ifsave = false;
verbose = true;

load(sprintf('%sdat_%s_%s_%.0fto%.0f',datadir,station,network,gc_lims(1),gc_lims(2)));
addpath('matguts')

%% cluster analysis on events
Z = linkage([eqar.gcarcs.*sind(eqar.seazs),eqar.gcarcs.*cosd(eqar.seazs)]);
T = cluster(Z,'cutoff',clust_cutoff,'criterion','distance');
nT = histcounts(T);
iclusts = find(nT>min4clust);

if verbose
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
end 

%% loop through clusters
for icl = 1:length(iclusts)

iclust = iclusts(icl);
gcarc_av = nan(length(RFphase),1);
seaz_av = nan(length(RFphase),1);
rayp_av = nan(length(RFphase),1);
norids = nan(length(RFphase),1);
evinds_save = cell(length(RFphase),1);

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
    evinds_clust = evinds;
    plot_rec_sec_align(eqar,RFphase{ip}(1),eqar.components{chp},filtfs,evinds)
%     clone_figure(66,67)
    %% Make matrix of data to cross-correlate        
    % grab the sample rate
    samprate = unique(round(1./diff(eqar.tt(:,:,ip))));
% %     % filter before cross correlation
% %     xcordat = eqar.dataZRT(:,:,chidx,ip);
% %     xcordat = filt_quick(xcordat,xcor_filtfs(1,ip),xcor_filtfs(2,ip),1./samprate,npoles);
% %     % taper extreme ends
% % %     xcordat = flat_hanning_win(xcortt,xcordat,xcortt(1)+5,xcortt(end)-5,5);
% %     % detrend
% %     xcordat = detrend(xcordat);
% %     % timing
% %     xcortt = eqar.tt(:,1,ip);
% %     xcinds = xcortt>=xcor_win(1) & xcortt<xcor_win(2); % indices of times within window
% %     % taper
% %     xcortt = eqar.tt(:,1,ip);
% %     xcordat = flat_hanning_win(xcortt,xcordat,xcortt(1)+5,xcortt(end)-5,5);
% %     
% %     % SNR after the filter
% %     nswin = find(xcortt>(max([datwind(1)+tapertime,-55])) & xcortt<=-5); % noise window 
% %     dwin  = find(xcortt>=-5 & xcortt<=10); % generous first arrival window 
% %     SNR_est = max(abs(xcordat(dwin,:)))./2./rms(xcordat(nswin,:));
% % 
% %     % window to xcorr time segment
% %     xcordat = detrend(xcordat(xcinds,:),'constant');
% %     
% %     % normalise
% %     for ie = 1:eqar.norids,
% %         dd = xcordat(:,ie);
% %         ddmax = max(abs(dd));
% %         xcordat(:,ie) = dd./ddmax; % flip and scale 
% % %         xcordat(:,ie) = dd.*ddsgn; % flip based on sign
% % %         xcordat(:,ie) = dd./maxab(dd); % flip and scale 
% %     end
% %     % window before cross correlation to xcor_window around arrival
% %     xcordat = flat_hanning_win(xcortt(xcinds),xcordat,xcor_win(1),xcor_win(2),1);
    %% SNRest
    % data clean for SNRest
    cp = struct('samprate',samprate,'fhi',xcor_filtfs(2,ip),'flo',xcor_filtfs(1,ip),...
                'pretime',-eqar.tt(1,1,ip),'prex',-datwind(1),'postx',datwind(2),...
                'taperx',1./diff(datwind),'npoles',npoles,'norm',1);
    [ SNRdat,~,~,~,~,SNRtt] = data_clean(  eqar.dataZRT(:,:,chidx,ip),cp ); 
    
    % SNR after the filter
    nswin = find(SNRtt>(max([datwind(1)+tapertime,-55])) & SNRtt<=-5); % noise window 
    dwin  = find(SNRtt>=-5 & SNRtt<=10); % generous first arrival window 
    SNR_est = max(abs(SNRdat(dwin,:)))./2./rms(SNRdat(nswin,:));
    % kill low SNRs
    evinds(SNR_est(evinds)<SNRmin) = [];
    plot_rec_sec_align(eqar,RFphase{ip}(1),eqar.components{chp},xcor_filtfs,evinds)

    %% polarity estimate
    % data clean for polest
    cp = struct('samprate',samprate,'fhi',10,'flo',1/100,...
                'pretime',-eqar.tt(1,1,ip),'prex',1,'postx',5,...
                'taperx',1./diff(xcor_win),'npoles',npoles,'norm',1);
    [ poldat,~,~,~,~,~] = data_clean(  eqar.dataZRT(:,:,chidx,ip),cp ); 
%     indtemplate = evinds(mindex(-SNR_est(evinds)));
%     t0template = STA_LTA(poldat(:,indtemplate),1./samprate,poltt);
%     iarr_templ = poltt<t0template+5;
    
    ddpol = zeros(eqar.norids,1);
    for ie = 1:eqar.norids
%         if ~ismember(ie,evinds), continue; end    
%         t0i = STA_LTA(poldat(:,ie),1./samprate,poltt,2.5);
%         if isnan(t0i), continue; end
%         iarr_i = poltt<t0i+2;
        
        % compare each indiv. to stack's polarity (as loose ref waveform)
        ddpol(ie) = polarity_est(sum(poldat,2),poldat(:,ie),5*samprate);
%         ddpol(ie) = polarity_est(poldat(iarr_templ,indtemplate),poldat(iarr_i,ie),5*samprate);
        % flip depending on polarity of maximum
%         xcordat(:,ie) = xcordat(:,ie).*ddpol(ie);
%     [ie,ddpol(ie)]

    end
    % kill arrs with no discernible pol
    evinds(ddpol(evinds)==0) = [];
    plot_rec_sec_align(eqar,RFphase{ip}(1),eqar.components{chp},xcor_filtfs,evinds,zeros(size(evinds)),ddpol)

     %% data clean for xcorr
    cp = struct('samprate',samprate,'fhi',xcor_filtfs(2,ip),'flo',xcor_filtfs(1,ip),...
                'pretime',-eqar.tt(1,1,ip),'prex',-xcor_win(1),'postx',xcor_win(2),...
                'taperx',1./diff(xcor_win),'npoles',npoles,'norm',1);
    [ xcordat,~,~,~,~,xcortt,~ ] = data_clean(  eqar.dataZRT(:,:,chidx,ip),cp ); 
     for ie = 1:eqar.norids,
        dd = xcordat(:,ie);
        ddmax = max(abs(dd));
        xcordat(:,ie) = dd.*ddpol(ie)./ddmax; % flip and scale 
%         xcordat(:,ie) = dd.*ddsgn; % flip based on sign
%         xcordat(:,ie) = dd./maxab(dd); % flip and scale 
     end
    
%% SNR after the filter
%     nswin = find(xcortt>(max([datwind(1)+tapertime,-55])) & xcortt<=-5); % noise window 
%     dwin  = find(xcortt>=-5 & xcortt<=10); % generous first arrival window 
%     SNR_est = max(abs(xcordat(dwin,:)))./2./rms(xcordat(nswin,:));


    %% CROSS CORRELATE - iteratively throw out low acor
    acor = -1;
%     evinds(SNR_est(evinds)<2) = [];
    if length(evinds)<3, continue; end
    
    if verbose
    figure(2),clf;plot(xcortt,xcordat(:,evinds))
    end
    
    while any(acor<0)
        fprintf('xcorr-ing %.0f traces\n',length(evinds))
        [dcor, dcstd, dvcstd, acor]=xcortimes(xcordat(:,evinds),1./samprate, -xcor_win(1), maxlag,verbose);
        fprintf('killing %.0f traces\n',sum(acor<0))
        evinds(acor<0) = [];
    end
    while any(acor<acormin)
        fprintf('xcorr-ing %.0f traces\n',length(evinds))
        
        [dcor, dcstd, dvcstd, acor]=xcortimes(xcordat(:,evinds),1./samprate, -xcor_win(1), maxlag,verbose);
        fprintf('killing %.0f traces\n',sum(acor<acormin))
        evinds(acor<acormin) = [];
    end
    
    if verbose
    plot_rec_sec_align(eqar,RFphase{ip}(1),eqar.components{chp},filtfs,evinds,-dcor,ddpol)
    end
    
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
        dataZRT_stk1(:,ic,ip) = sum(dataZRT_proc(:,evinds,ic,ip),2)/length(evinds); %#ok<SAGROW>
    end
    

    %% Xcorr again!
    fprintf('Xcorr on stack\n')
    dataZRT_proc2 = nan([size(eqar.dataZRT(:,:,1,ip)),3]);
    acor = zeros(length(evinds_clust),ip);
    dcor = zeros(length(evinds_clust),ip);
    SNR_all2 = zeros(length(evinds_clust),3);
    for ie = 1:length(evinds_clust)
        evind = evinds_clust(ie);        
        dd_s = squeeze(eqar.dataZRT(:,evind,:,ip));
        % fix
%         dd_s = interp1(eqar.tt(:,evind,ip), dd ,eqar.tt(:,evind,ip)+dcor(ie)); % interp to new time axis with the shift
        dd_s(isnan(dd_s)) = 0;
        % filt
        dd_sf = filt_quick(dd_s,filtfs(1),filtfs(2),1./samprate);
        % taper
        dd_sft = flat_hanning_win(eqar.tt(:,evind,ip),dd_sf,datwind(1),datwind(2),tapertime);
        % normalise
        normval = max(max(abs(dd_sft)));
        % polarity
        [~,pol] = polarity_est(dataZRT_stk1(:,chp,ip),...
                           flat_hanning_win(eqar.tt(:,evind,ip),dd_sft(:,chp),-20,20,tapertime),...
                            5*samprate); % using polest_cc
        if pol==0, acor(ie)=0; continue; end
        dd_sftnp = dd_sft./normval./pol;
        % align
        [ddcor,~,~,aacor] = xcortimes([dd_sftnp(:,chp),dataZRT_stk1(:,chp,ip)],1./samprate,0,5,0);
        acor(ie,ip) = sum(aacor)-1;
        dcor(ie,ip) = -diff(ddcor);
        dataZRT_proc2(:,evind,:,ip) = interp1(eqar.tt(:,evind,ip), dd_sftnp ,eqar.tt(:,evind,ip)+dcor(ie,ip));

        %% SNR
        nswin = find(eqar.tt(:,evind,ip)<=-5 & eqar.tt(:,evind,ip)>(max([datwind(1)+tapertime,-55]))); % noise window 
        dwin  = find(eqar.tt(:,evind,ip)>=-2 & eqar.tt(:,evind,ip)<=10); % generous first arrival window 
        for ic = 1:3
            SNR_all2(ie,ic) = max(abs(detrend(dd_sftnp(dwin,ic))))./2./rms(detrend(dd_sftnp(nswin,ic)));
        end
%         figure(33), clf, set(gcf,'pos',[84 675 1156 423]);
%         subplot(221);plot(eqar.tt(:,evind,ip),dd_s(:,chp)/normval,eqar.tt(:,evind,ip),dataZRT_stk1(:,chp,ip))
%         subplot(222);plot(eqar.tt(:,evind,ip),dd_s(:,chd)/normval,eqar.tt(:,evind,ip),dataZRT_stk1(:,chd,ip))
%         subplot(223);plot(eqar.tt(:,evind,ip),dataZRT_proc2(:,evind,chp,ip),eqar.tt(:,evind,ip),dataZRT_stk1(:,chp,ip))
%         subplot(224);plot(eqar.tt(:,evind,ip),dataZRT_proc2(:,evind,chd,ip),eqar.tt(:,evind,ip),dataZRT_stk1(:,chd,ip))
%         pause
    end
    
    % get indices of events with good acor w/ the stack
%     gdinds = evinds_clust(acor(:,ip)>.7);
%     gdinds = evinds_clust(SNR_all2(:,chp)>SNRmin & SNR_all2(:,chd)>SNRmin/10);
    gdinds = evinds_clust(acor(:,ip)>acormin & SNR_all2(:,chp)>SNRmin & SNR_all2(:,chd)>SNRmin/10);
   	for ic = 1:3
        dataZRT_stk(:,ic,ip) = sum(dataZRT_proc2(:,gdinds,ic,ip),2)/length(gdinds); %#ok<SAGROW>
    end

    %% estimate Vp, Vs at surface from trace
    err_surface = [];
    for ie = 1:length(evinds)
    ei = evinds(ie);
    [ ~,~,err_surface(:,:,ie),vp_range,vs_range ] = ...
            surf_vp_vs_est( squeeze(dataZRT_proc(:,ei,:,ip))*[1 0 0;0 -1 0;0 0 1],eqar.tt(:,ei),eqar.rayps(ei),[0 10] );
    end   
    figure(1),clf, contourf(vs_range,vp_range,sum(err_surface,3),40,'edgecolor','none'), colorbar
    Esurf(:,:,ip) = sum(err_surface,3);
    [~,x,y] = mingrid(Esurf(:,:,ip));
    Vp_est(icl,ip) = vp_range(y);
    Vs_est(icl,ip) = vs_range(x);
   


    if verbose
    figure(34)
    subplot(211);plot(eqar.tt(:,1,ip),dataZRT_stk1(:,chp,ip),eqar.tt(:,1,ip),dataZRT_stk(:,chp,ip))
    subplot(212);plot(eqar.tt(:,1,ip),dataZRT_stk1(:,chd,ip),eqar.tt(:,1,ip),dataZRT_stk(:,chd,ip))
    end
 
	%% continue if not enough for stack
    if length(gdinds) < min4stack, fprintf('Not enough good traces\n'), continue; end
    fprintf('%.0f good traces\n',length(gdinds))
    
    %% automatic pick phase onset
    figure(20), clf, set(gcf,'pos',[15 225 1426 420]), hold on
    htr = plot(eqar.tt(:,1,ip),dataZRT_stk(:,:,ip),'linewidth',1.5);
    title('PICK PHASE ONSET','fontsize',22)
    set(gca,'xlim',[-50 50],'ylim',[-1. 1.],'fontsize',18)
    legend(eqar.components)
    
    % first use sta_lta to get approx phase onset
    istalta = eqar.tt(:,1,ip)>=-50 & eqar.tt(:,1,ip)<50;
    [stalta_t0] = STA_LTA( dataZRT_stk(istalta,chp,ip),1./samprate,eqar.tt(istalta,1,ip));
    % now, within 5s of sta_lta onset, find the max amp.
    imax = eqar.tt(:,1,ip)>=stalta_t0-5 & eqar.tt(:,1,ip)<stalta_t0+5;
    [~,tshft] = crossing(dataZRT_stk(imax,chp,ip),eqar.tt(imax,1,ip),maxab(dataZRT_stk(imax,chp,ip)));
    
    tt_stk(1:size(eqar.tt,1),ip) = eqar.tt(:,1,ip) - tshft(1);
    tt_stk = round_level(tt_stk,1/samprate);
    for ic = 1:length(htr)
        htr(ic).XData = tt_stk(:,ip);
        htr(ic).YData = dataZRT_stk(:,ic,ip);
    end

    pause(0.1)

    %% other averaged parms
    gcarc_av(ip) = mean(eqar.gcarcs(gdinds));
    seaz_av(ip) = mean(eqar.seazs(gdinds));
    rayp_av(ip) = mean(eqar.rayps(gdinds,ip));
    norids(ip,1) = length(gdinds);
    evinds_save{ip,1} = gdinds;
 
    %% nice plot
    figure(21), clf, set(gcf,'pos',[92 404 802 642])
    for ic = 1:length(eqar.components)
        subplot(3,1,ic), hold on
        plot(tt_stk(:,ip),dataZRT_proc(:,gdinds,ic,ip),'color',[.8 .89 .9],'linewidth',0.5)
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
    
% [ Vp_est_stk_p,Vs_est_stk_p,errmap_p ] = ...
%             surf_vp_vs_est( dataZRT_stk(:,:,1)*[1 0 0;0 -1 0;0 0 1],tt_stk(:,1),rayp_av(1),[-5 10],1 );
% [ Vp_est_stk_s,Vs_est_stk_s,errmap_s,vp_range,vs_range ] = ...
%             surf_vp_vs_est( dataZRT_stk(:,:,2)*[1 0 0;0 -1 0;0 0 1],tt_stk(:,2),rayp_av(2),[-5 10],1 );
% errmap = errmap_p/mingrid(errmap_p) + errmap_s/mingrid(errmap_s);
% [~,x,y] = mingrid(errmap);
% Vp_est_stk = vp_range(y)
% Vs_est_stk = vs_range(x)
% figure(5),clf, hold on; contourf(vs_range,vp_range,errmap,40,'edgecolor','none'), colorbar

if ~exist('rayp_av','var'), continue; end

%% kill ones with too few evts in each stack
keepstk = norids >= min4stack;
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
              'dataZRT',dataZRT_stk,'tt',tt_stk,'stackmethod','doublestack');
          
if ifsave
    fprintf('SAVING\n')
    arfile = sprintf('avar_dat_%s_%s_%.0f_%.0f',station,network,mean(gcarc_av),mean(seaz_av));
    save(arfile,'avar');
    if exist('AVARS','dir'); movefile([arfile,'.mat'],'AVARS'); end
else
    pause
end

end % loop on clusters
    

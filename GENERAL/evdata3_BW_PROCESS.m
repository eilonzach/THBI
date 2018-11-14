close all
clear 
%% Setup
proj = struct('name','NWUS');
proj.dir = ['~/Documents/MATLAB/BayesianJointInv/',proj.name];
wd = pwd;
addpath(wd);

%% load project, station, and request details and request details
load([proj.dir,'/project_details.mat']);
load([proj.infodir,'stations.mat']);
load([proj.infodir,'/data_request_details.mat']);

generation = 30; % generation of solution and data processing

phases = {'Ps','Sp'}; % {'Ps','Sp'}

ifsave = true;
ifcalcsurfV = true; % option to re-calculate surface 
ifPSV = true;
ifverbose = false;
ifpatricide = true;

% FINAL DATA processing parms
%=================================
gc_lims = [30,75]; % 
filtfs = [1/50;1000]; % filter frequencies pre PCA, in Hz
% dat_win = [-70 70]; % in sec around main arrival
Npass = 1; % 1 (causal) 2 (acausal)
%=================================

% PATRICIDE - kill parts of waveforms we think are incorrectly-transformed parent...
% P-dat
patricide_win(:,:,1) = [-3   35  % P  comp
                     0   35];% SV comp     % was [1 35]
% S-dat
patricide_win(:,:,2) = [-35   0  % P  comp  % was [-35 -1]
                     -35   5];% SV comp   

% polest parms
pol_filtfs = [1/40 1];
pol_win = [-5 5];

% Processing parms
npoles = 2;
tapertime = 5; % s at beginnning and end of window to taper over

% SNR parms
SNR_filtfs = [0.02;100];
% SNR_filtfs = [1/100;1000];
SNR_win = [[-85 -15];[-10 10]]; % noise window, data window, in sec from main arrival
SNRmin = 2;%3!!

% P/SV parms
psv_SNRmin = 3;%3!!
psv_SNRmax = 30;% for max weighting of error surfaces
% psv_filtfs = [[0.01;1],[0.066667;0.5]]; % p,s
psv_filtfs = [[0.4;2],[0.05;1/1]]; % p,s
psv_win = [[-10;10],[-10;10]]; % p,s  [[-5;5],[-3;6]]

Vp_default = 3.1*1.8;%RSSD: 5.0; % WVOR: 3.45 % HYB: 4.5
Vs_default = 3.1;%RSSD: 3.0; % WVOR: 2.46 % HYB: 3.18

% Xcorr parms
xcor_filtfs = [[0.4;2],[0.05;1/1]]; % p,s
xcor_win = [-20 5];
acormin = 0.60; % 0.7
stkacormin = 0.80; % 0.9
maxlag = 10; % in seconds

% cluster parms
maxdepth = 200;  % maximum EQ depth
min4clust = 25; % minimum # of EQ to even try the stack
clust_bazwin = 30;
clust_arcwin = 20;
clust_cutoff = 5; % minimum Delta of clustered EQ locations (deg)
min4stack = 10;  % minumum # of traces in order to keep the stack

% addpath('~/Documents/MATLAB/BayesianJointInv/US/matguts/')
addpath('~/Documents/MATLAB/BayesianJointInv/functions/')

% load Shen & Ritzwoller model for migration
ncfile = '~/Work/data/models_seismic/SR16_3d_Vs/US.2016.nc';
[ sr16_model ] = read_seismodel_nc( ncfile );


%% ==================  LOOP OVER STATIONS IN DB  ================== 
%% ==================  LOOP OVER STATIONS IN DB  ================== 
for is = 28:stainfo.nstas
station = stainfo.stas{is};
network = stainfo.nwk{is};
fprintf( '\n --------- STATION %s,  NETWORK %s  --------- \n',station,network);

rawfile = dir(sprintf('%sdat_%s_%s*',proj.rawdatadir,station,network));
if isempty(rawfile)
    fprintf( '>>> No data for sta=%s, nwk=%s <<<\n',station,network);
    continue
end
if length(rawfile)>1, error('More than one data file for this sta+nwk...'); end


%% Load data, prep outputs
load([rawfile.folder,'/',rawfile.name]);
odir = [proj.STAinversions,station,'_',network,'_dat',num2str(generation)]; 
if ~exist(odir,'dir')~=7, mkdir(odir); end

%% get Shen and Ritzwoller model for migration; tack ak135 onto the bottom of it
vs = squeeze(sr16_model.Vsv(mindex(sr16_model.lon,mod(stainfo.slons(is),360)),mindex(sr16_model.lat,stainfo.slats(is)),:));
vp = squeeze(sr16_model.Vpv(mindex(sr16_model.lon,mod(stainfo.slons(is),360)),mindex(sr16_model.lat,stainfo.slats(is)),:));
Z = sr16_model.Z;
if mean(vs)>10
    addpath('~/Documents/MATLAB/seizmo/models/')
    akmod = ak135('depths',[0:5:300]');
    Z = akmod.depth; 
    vs = akmod.vs; 
    vp = akmod.vp; 
elseif max(Z)<300 
    addpath('~/Documents/MATLAB/seizmo/models/')
    akmod = ak135('depths',[max(Z)+5:5:300]');
    Z = [Z;akmod.depth]; 
    vs = [vs;akmod.vs]; 
    vp = [vp;akmod.vp]; 
end % extend all models to 400 km depth
Vmodel = struct('z',Z,'VP',vp,'VS',vs);

%% Estimate Vp and Vs in the crust by minimising P/SV
% Use all data available
% weight error maps by SNR
if ifcalcsurfV && ifPSV
    [Vp_surf,Vs_surf] = evdata_VpVs_est(eqar,phases,psv_filtfs,psv_win,[],[psv_SNRmin,psv_SNRmax],SNR_filtfs,SNR_win,ifverbose);
    Vp_surf = Vp_surf.comp3_wtstk;
    Vs_surf = Vs_surf.comp3_wtstk;
    pause(0.1)
else
    Vp_surf = Vp_default;
    Vs_surf = Vs_default;
end

SNR_est = zeros(eqar.norids,2,2); % orid,component,phase

% return

%% Cluster analysis over back-azimuth
% Z = linkage([eqar.gcarcs.*sind(eqar.seazs),eqar.gcarcs.*cosd(eqar.seazs)]);
% T = cluster(Z,'cutoff',clust_cutoff,'criterion','distance');
% nT = histcounts(T);
% iclusts = find(nT>=min4clust);

[ iclusts,nclust ] = eqcluster( eqar.seazs,eqar.gcarcs,min4clust,clust_bazwin,clust_arcwin,1 );
T = iclusts;
iclusts = 1:nclust;

% plot_clust_radplot( eqar.seazs,eqar.gcarcs,T )

%% Subset to cluster with largest N (can relax this later to loop on all)
for icl = 1:length(iclusts)
    fprintf('\n==================================================\n')
    fprintf('===== Cluster %.0f  gcarc ~ %4.1f  seaz ~ %5.1f  =====\n\n',icl,mean(eqar.gcarcs(T==icl)),mean(eqar.seazs(T==icl)))
    
    plot_dat = cell({});
    nevinds = zeros(length(phases),1);
    evinds_use = cell(length(phases),1);
    rayp_av = nan(length(phases),1);
    gcarc_av = nan(length(phases),1);
    seaz_av = nan(length(phases),1);
    data_stk =[];
    tt_stk = [];

% Individual surface transform velocities...    
% if ifcalcsurfV && ifPSV
%     evinds_do_psv  = [1:eqar.norids]';
%     evinds_do_psv(~(T==iclusts(icl))) = [];    
% 	[Vp_surf,Vs_surf] = evdata_VpVs_est(eqar,phases,psv_filtfs,psv_win,evinds_do_psv,[psv_SNRmin,psv_SNRmax],SNR_filtfs,SNR_win,verbose);
%     Vp_surf = Vp_surf.comp3_wtstk;
%     Vs_surf = Vs_surf.comp3_wtstk;
% else
%     Vp_surf = Vp_default;
%     Vs_surf = Vs_default;
% end


    
%% loop on phases   
for ip = 1:length(phases)
    fprintf('Processing %s data, cluster %.0f\n',phases{ip},icl)
    evinds  = [1:eqar.norids]';
    evinds(~(T==iclusts(icl))) = [];
    evinds(eqar.gcarcs(evinds)>max(gc_lims) | eqar.gcarcs(evinds)<min(gc_lims)) = [] ;
    evinds(eqar.edeps(evinds)>maxdepth) = [] ;
    if strcmp(phases{ip},'Sp'),evinds(eqar.gcarcs(evinds)<60) = []; end
    if isempty(evinds) || length(evinds)<min4stack, continue; end
    
    %% QC step to kill anything with gradients that are way way out of the ordiary (e.g. from instrument spike)
    evinds(any(abs(diff(eqar.dataZRT(:,evinds,1,ip)))>30*mean(rms(diff(eqar.dataZRT(:,evinds,1,ip))))))=[];
    evinds(any(abs(diff(eqar.dataZRT(:,evinds,2,ip)))>30*mean(rms(diff(eqar.dataZRT(:,evinds,2,ip))))))=[];
    evinds(any(abs(diff(eqar.dataZRT(:,evinds,1,ip)))>50*mean(rms(diff(eqar.dataZRT(:,evinds,1,ip))))))=[];
    evinds(any(abs(diff(eqar.dataZRT(:,evinds,2,ip)))>50*mean(rms(diff(eqar.dataZRT(:,evinds,2,ip))))))=[];
    fprintf(' > %.0f evinds\n',length(evinds))
    
    %% grab the sample rate
    samprate = unique(round(1./diff(eqar.tt(:,:,ip))));
    ph = phases{ip};
    switch ph(1)
        case 'P', chp = 1; chd = 2; tdir = +1;% P parent, SV daughter
        case 'S', chp = 2; chd = 1; tdir = -1;% SV parent, P daughter
    end    
    
    % ----------------------- save data level 1 -----------------------
    plot_dat{1,ip} = struct('data',eqar.dataZRT(:,evinds,1:2,ip),'tt',eqar.tt(:,1:2,ip));
    % ----------------------- save data level 1 -----------------------
        
%% Transform to P-SV
    if ifPSV, fprintf(' > transform to P-SV coord system\n'), end
    datPSV = zeros(size(eqar.dataZRT,1),eqar.norids,2);
    
    for ie = 1:length(evinds)
        evind = evinds(ie);
        if ifPSV
        % Z already positive down...
        [datP,datSV] = ...
            Rotate_XZ_to_PSV(eqar.dataZRT(:,evind,2,ip),eqar.dataZRT(:,evind,1,ip),...
                             Vp_surf,Vs_surf,rayp_sdeg2skm(eqar.rayps(evind,ip)));
        datPSV(:,evind,1) = datP ./maxgrid(abs([datP,datSV]));
        datPSV(:,evind,2) = datSV./maxgrid(abs([datP,datSV]));
        else
        datPSV(:,evind,1) = -eqar.dataZRT(:,evind,1,ip)./maxgrid(abs(eqar.dataZRT(:,evind,1:2,ip)));
        datPSV(:,evind,2) =  eqar.dataZRT(:,evind,2,ip)./maxgrid(abs(eqar.dataZRT(:,evind,1:2,ip)));
        end            
    end
    
    % high pass
    for ie = 1:length(evinds)
        evind = evinds(ie);
        switch ip
            case 1
                datPSV(:,evind,1) = filt_quick(datPSV(:,evind,1),filtfs(1),filtfs(2),1/samprate,npoles,Npass,-1);
                datPSV(:,evind,2) = filt_quick(datPSV(:,evind,2),filtfs(1),filtfs(2),1/samprate,npoles,Npass,-1);
            case 2
                datPSV(:,evind,1) = filt_quick(datPSV(:,evind,1),filtfs(1),filtfs(2),1/samprate,npoles,Npass);
                datPSV(:,evind,2) = filt_quick(datPSV(:,evind,2),filtfs(1),filtfs(2),1/samprate,npoles,Npass);
        end
    end


    % ----------------------- save data level 2 -----------------------
    plot_dat{2,ip} = struct('data',datPSV(:,evinds,1:2),'tt',eqar.tt(:,evinds,ip));
    % ----------------------- save data level 2 -----------------------
    

%% Calculate SNR for all arrivals, for both components
    fprintf(' > estimating SNR\n')
    % data clean for SNRest
    cp_snr = struct('samprate',samprate,'fhi',SNR_filtfs(2),'flo',SNR_filtfs(1),...
                'pretime',-eqar.tt(1,1,ip),'prex',tapertime-SNR_win(1,1),'postx',SNR_win(2,2)+tapertime,...
                'npoles',npoles,'norm',1);
    cp_snr.taperx = tapertime./(cp_snr.prex+cp_snr.postx);
    [ SNRdatp,~,~,~,~,SNRtt] = data_clean(  datPSV(:,evinds,chp),cp_snr ); 
    [ SNRdatd,~,~,~,~,SNRtt] = data_clean(  datPSV(:,evinds,chd),cp_snr ); 

    % SNR after the filter
    nswin = find(SNRtt>=SNR_win(1,1) & SNRtt < SNR_win(1,2)); % noise window 
    dwin  = find(SNRtt>=SNR_win(2,1) & SNRtt < SNR_win(2,2)); % generous first arrival window 
    SNR_est(evinds,chp,ip) = max(abs(SNRdatp(dwin,:)))./2./rms(detrend(SNRdatp(nswin,:)));
    SNR_est(evinds,chd,ip) = max(abs(SNRdatd(dwin,:)))./2./rms(detrend(SNRdatd(nswin,:)));
        
%% QC on SNR
    fprintf(' > eliminate %.0f low SNR\n',sum(SNR_est(evinds,chp,ip)<SNRmin))
    evinds(SNR_est(evinds,chp,ip)<SNRmin) = [];
    if length(evinds)<min4stack, continue; end


%% Polarity estimate
    fprintf(' > correcting polarity\n')
    % data clean for polest
    cp = struct('samprate',samprate,'fhi',pol_filtfs(2),'flo',pol_filtfs(1),...
                'pretime',-eqar.tt(1,1,ip),'prex',-pol_win(1),'postx',pol_win(2),...
                'taperx',0.5./diff(pol_win),'npoles',npoles,'norm',1);
    [ poldat,~,~,~,~,~] = data_clean(  datPSV(:,evinds,chp),cp ); 
    
    % make preliminary stack to compare polarity
    [ B,ind ] = maxab( poldat ); % find maximum points 
    ind = ind - min(ind);
    poldat_shiftscale = nan(size(poldat,1)+max(ind),size(poldat,2));
    for ii = 1:length(evinds)
    % shift and scale all traces to align these max points (flipping
    % reverse polarity ones just for this prelim stack - not in data, yet)
        poldat_shiftscale(:,ii) = [zeros(max(ind)-ind(ii),1);poldat(:,ii)./B(ii);zeros(ind(ii),1)];
    end
    prelimstack = sum(poldat_shiftscale(floor(max(ind)/2):end-1-ceil(max(ind)/2),:),2); % cut stack to be same length as data
    
    ddpol = zeros(eqar.norids,1);
    for ie = 1:length(evinds)
        evind = evinds(ie);
        % compare each indiv. to stack's polarity (as loose ref waveform)
        ddpol(evind) = polarity_est(prelimstack,poldat(:,ie),5*samprate);
    end
    % kill arrs with no discernible pol
    fprintf('      killing %.0f traces with bad pol\n',sum(ddpol(evinds)==0))
    evinds(ddpol(evinds)==0) = [];
    if length(evinds)<min4stack, continue; end

    for ie = 1:length(evinds)
        evind = evinds(ie);
        datPSV(:,evind,:) = datPSV(:,evind,:)*ddpol(evind);
    end
    
    % ----------------------- save data level 3 -----------------------
    plot_dat{3,ip} = struct('data',datPSV(:,evinds,1:2),'tt',eqar.tt(:,evinds,ip));
    % ----------------------- save data level 3 -----------------------


%% Cross correlate (inc. filter to xcorr band)
    fprintf(' > Xcorr parent\n')
    cp_xcor = struct('samprate',samprate,'fhi',xcor_filtfs(2,ip),'flo',xcor_filtfs(1,ip),...
                'pretime',-eqar.tt(1,1,ip),'prex',-xcor_win(1),'postx',xcor_win(2),...
                'taperx',1./diff(xcor_win),'npoles',npoles,'norm',1);
    [ xcordat,~,~,~,~,xcortt] = data_clean(  datPSV(:,evinds,chp),cp_xcor ); 

    % CROSS CORRELATE - iteratively throw out low acor
    acor = -1;
    % get rid of negative acors
    while any(acor<0)
        fprintf('    xcorr-ing %.0f traces\n',size(xcordat,2))
        [dcor, dcstd, dvcstd, acor]=xcortimes(xcordat,1./samprate, -xcor_win(1), maxlag,0);
        fprintf('      killing %.0f traces\n',sum(acor<0))
        xcordat(:,acor<0) = [];
        evinds(acor<0) = [];
    end
    % get rid of acors below acormin threshold
    while any(acor<acormin) && length(evinds)>=min4stack
        fprintf('    xcorr-ing %.0f traces\n',size(xcordat,2))
        [dcor, dcstd, dvcstd, acor]=xcortimes(xcordat,1./samprate, -xcor_win(1), maxlag,0);
        fprintf('      killing %.0f traces\n',sum(acor<acormin))
        xcordat(:,acor<acormin) = [];
        evinds(acor<acormin) = [];
    end
    if length(evinds)<min4stack, continue; end
    
    %% get average time of maximum 
    Tmax = zeros(length(evinds),1);
    for ie = 1:length(evinds)
        evind = evinds(ie);
        ind = abs(eqar.tt(:,evind,ip))<6;
        [~,Tmax(ie)] = crossing(datPSV(ind,evind,chp),eqar.tt(ind,evind,ip)+dcor(ie),max(datPSV(ind,evind,chp)));
    end
    %use this to correct the dcors give max at origin
    dcor = dcor+mean(Tmax);
    
    
    %% Align traces, migrate, and stack
    fprintf(' > align traces\n')
    CdatPSV = nan(size(datPSV,1),eqar.norids,2);
    for ie = 1:length(evinds)
        evind = evinds(ie);
        % align
        CdatPSV(:,evind,1) = interp1(eqar.tt(:,evind,ip), datPSV(:,evind,1) ,eqar.tt(:,evind,ip)+dcor(ie),'linear',0); % interp to new time axis with the shift
        CdatPSV(:,evind,2) = interp1(eqar.tt(:,evind,ip), datPSV(:,evind,2) ,eqar.tt(:,evind,ip)+dcor(ie),'linear',0); % interp to new time axis with the shift
        % 100s high pass
%         CdatPSV(:,evind,1) = filt_quick(CdatPSV(:,evind,1),filtfs(1),filtfs(2),1/samprate,npoles,Npass);
%         CdatPSV(:,evind,2) = filt_quick(CdatPSV(:,evind,2),filtfs(1),filtfs(2),1/samprate,npoles,Npass);
    end

    %% migrate
    % get ray parameter of cluster
    raypavs_clust = rayp_sdeg2skm(mean(eqar.rayps(evinds,ip)));
    CdatPSVmig = CdatPSV;
    for ie = 1:length(evinds)
        evind = evinds(ie);
        posTdir = -1i^(2*ip);% 1 for ip==1, -1 for ip==2;
        indmig = find(posTdir*eqar.tt(:,evind,ip) >= 0);
        indpre = find(posTdir*eqar.tt(:,evind,ip) < 0);
        
        ott = migrate_PorS_conv( eqar.tt(indmig,evind,ip),Vmodel,rayp_sdeg2skm(eqar.rayps(evind,ip)),raypavs_clust,phases{ip});
        
        indmig = indmig(~isnan(ott));
        % resolve onto new time axis
        migdat1 = interp1(ott(~isnan(ott)),CdatPSV(indmig,evind,1),eqar.tt(:,evind,ip),'linear',nan);
        migdat2 = interp1(ott(~isnan(ott)),CdatPSV(indmig,evind,2),eqar.tt(:,evind,ip),'linear',nan);
        % add in pre-arrival data
        migdat1(indpre) = CdatPSV(indpre,evind,1);
        migdat2(indpre) = CdatPSV(indpre,evind,2);
        
        CdatPSVmig(:,evind,1) = migdat1;
        CdatPSVmig(:,evind,2) = migdat2;
    end
        
    
%     figure;plot(eqar.tt(:,1,ip),CdatPSVmig(:,evinds,1))
%     figure;plot(eqar.tt(:,1,ip),CdatPSV(:,evinds,1))
    
    % STACK #1
    fprintf(' > stack #1\n')
    SdatPSV = squeeze(sum(CdatPSVmig(:,evinds,:),2));    

    %     figure;plot(eqar.tt(:,1,ip),SdatPSV,eqar.tt(:,1,ip),SdatPSVmig)
    
    % ----------------------- save data level 4 -----------------------
    plot_dat{4,ip} = struct('data',CdatPSVmig(:,evinds,1:2),'tt',eqar.tt(:,evinds,ip));
    % ----------------------- save data level 4 -----------------------

    % ----------------------- save data level 5 -----------------------
    plot_dat{5,ip} = struct('data',SdatPSV(:,1:2),'tt',eqar.tt(:,1,ip));
    % ----------------------- save data level 5 -----------------------
          
        
    %% Xcor indiv traces with stack
    fprintf(' > Xcorr on stack\n')
    [ xcorSTK] = data_clean(SdatPSV(:,chp),cp_xcor ); 
    acor_stk = zeros(length(evinds),1); dcor_stk = zeros(length(evinds),1);
    for ie = 1:length(evinds)
        [dd,~,~,aa] = xcortimes([xcorSTK,xcordat(:,ie)],1./samprate,-xcor_win(1),maxlag,0);
        acor_stk(ie) = sum(aa)-1;
        dcor_stk(ie) = diff(dd);
    end
    % reject indiv traces that don't look enough like the stack...
    fprintf('      killing %.0f traces that do not look like stack\n',sum(acor_stk<stkacormin))
	evinds(acor_stk<stkacormin) = [];
    dcor_stk(acor_stk<stkacormin) = [];
    if length(evinds)<min4stack, continue; end
    
    %% Re-align traces and re-stack
    fprintf(' > Re-align traces\n')
    C2datPSV = nan(size(datPSV,1),eqar.norids,2);
    for ie = 1:length(evinds)
        evind = evinds(ie);
        % align
        C2datPSV(:,evind,1) = interp1(eqar.tt(:,evind,ip), datPSV(:,evind,1) ,eqar.tt(:,evind,ip)+dcor_stk(ie),'linear',0); % interp to new time axis with the shift
        C2datPSV(:,evind,2) = interp1(eqar.tt(:,evind,ip), datPSV(:,evind,2) ,eqar.tt(:,evind,ip)+dcor_stk(ie),'linear',0); % interp to new time axis with the shift
    end
    
    %% migrate2
    % get ray parameter of cluster
    raypavs_clust = rayp_sdeg2skm(mean(eqar.rayps(evinds,ip)));
    C2datPSVmig = C2datPSV;
    for ie = 1:length(evinds)
        evind = evinds(ie);
        posTdir = -1i^(2*ip);% 1 for ip==1, -1 for ip==2;
        indmig = find(posTdir*eqar.tt(:,evind,ip) >= 0);
        indpre = find(posTdir*eqar.tt(:,evind,ip) < 0);
        
        ott = migrate_PorS_conv( eqar.tt(indmig,evind,ip),Vmodel,rayp_sdeg2skm(eqar.rayps(evind,ip)),raypavs_clust,phases{ip});
        
        indmig = indmig(~isnan(ott));
        % resolve onto new time axis
        migdat1 = interp1(ott(~isnan(ott)),C2datPSV(indmig,evind,1),eqar.tt(:,evind,ip),'linear',nan);
        migdat2 = interp1(ott(~isnan(ott)),C2datPSV(indmig,evind,2),eqar.tt(:,evind,ip),'linear',nan);
        % add in pre-arrival data
        migdat1(indpre) = C2datPSV(indpre,evind,1);
        migdat2(indpre) = C2datPSV(indpre,evind,2);
        
        C2datPSVmig(:,evind,1) = migdat1;
        C2datPSVmig(:,evind,2) = migdat2;
    end
    
    % STACK #2
    fprintf(' > Re-stack\n')
    S2datPSV = squeeze(sum(C2datPSVmig(:,evinds,:),2));
    S2datPSV = S2datPSV./maxab(S2datPSV(:));
    Stt = eqar.tt(:,1,ip);

    % ----------------------- save data level 6 -----------------------
    plot_dat{6,ip} = struct('data',C2datPSVmig(:,evinds,1:2),'tt',eqar.tt(:,evinds,ip));
    % ----------------------- save data level 6 -----------------------

    % ----------------------- save data level 7 -----------------------
    plot_dat{7,ip} = struct('data',S2datPSV(:,1:2),'tt',Stt);
    % ----------------------- save data level 7 -----------------------

    %% automatic pick phase onset
    fprintf(' > t0 main phase\n')
    figure(20+10*ip), clf, set(gcf,'pos',[15 (220+300*(ip-1)) 1426 420]), hold on
    htr = plot(Stt,S2datPSV,'linewidth',1.5);
    set(gca,'xlim',[-50 50],'ylim',[-1. 1.],'fontsize',18)
    legend({'P','SV'})
    
%     [~,tshftu] = crossing(S2datPSV(:,chp),Stt,max(S2datPSV(:,chp)));
%     [~,tshftd] = crossing(-S2datPSV(:,chp),Stt,max(-S2datPSV(:,chp)));
%     tshft = min([tshftd,tshftu]);
    [~,tshft] = crossing(S2datPSV(:,chp),Stt,max(S2datPSV(:,chp)));
    title(sprintf('Cluster %.0f  Phase %s $\\Delta\\sim%4.1f$  $\\phi\\sim%5.1f$ \n\n',icl,phases{ip},mean(eqar.gcarcs(evinds)),mean(eqar.seazs(evinds))),'fontsize',22,'interpreter','latex')

    
    Stt = Stt - tshft;
    Stt = round_level(Stt,1/samprate);
    for ic = 1:length(htr)
        htr(ic).XData = Stt;
        htr(ic).YData = S2datPSV(:,ic);
    end

    pause(0.1)
    if ifverbose
        fprintf('>> Click/hit space to continue <<\n');
        pause
    end
    % ----------------------- save data level 8 -----------------------
    plot_dat{8,ip} = struct('data',S2datPSV(:,1:2),'tt',Stt);
    % ----------------------- save data level 8 -----------------------
    
    %% PATRICIDE by killing parts of waveforms we think are un-rotated parent...
    if ifpatricide

    for ic = 1:2
        S2datPSV_patricide(:,ic) = flat_hanning_win(Stt,S2datPSV(:,ic),patricide_win(ic,1,ip),patricide_win(ic,2,ip),2);
    end
    
    figure(21+10*ip), clf, set(gcf,'pos',[15 425 1426 420]), hold on
    htr = plot(Stt,S2datPSV_patricide,'linewidth',1.5);
    set(gca,'xlim',[-50 50],'ylim',[-1. 1.],'fontsize',18)
    legend({'P','SV'})
    
    % sub in patricidal data
    S2datPSV = S2datPSV_patricide;
    end
    
    %% SNR on each component
    fprintf(' > Estimate stack SNR\n')
    nswin = find(Stt>=SNR_win(1,1) & Stt < SNR_win(1,2)); % noise window 
    dwin  = find(Stt>=2*SNR_win(2,1) & Stt < 2*SNR_win(2,2)); % generous first arrival window 
    SNR_stk(chp,ip) = max(abs(S2datPSV(dwin,chp)))./2./rms(detrend(S2datPSV(nswin,chp)));
    SNR_stk(chd,ip) = max(abs(S2datPSV(dwin,chd)))./2./rms(detrend(S2datPSV(nswin,chd)));
    SNR_stk(3,ip) = 0;

    %% Save prep
    if length(evinds)>=min4stack
    nevinds(ip) = length(evinds);
    evinds_use{ip} = evinds;
    rayp_av(ip) = mean(eqar.rayps(evinds,ip));
    gcarc_av(ip) = r2d(mean_ang(d2r(eqar.gcarcs(evinds))));
    seaz_av(ip) = r2d(mean_ang(d2r(eqar.seazs(evinds))));
    data_stk(:,:,ip) = [S2datPSV,zeros(length(Stt),1)];
    tt_stk(:,ip) = Stt;
    else
    nevinds(ip) = 0;
    evinds_use{ip} = 0;
    rayp_av(ip) = nan;
    gcarc_av(ip) = nan;
    seaz_av(ip) = nan;
    data_stk(:,:,ip) = nan;
    tt_stk(:,ip) = nan;
    end
    
end % on phase loop
if isempty(plot_dat) || size(plot_dat,1)<8
    continue
end

% TRY TO BACK OUT SURFACE VELOCITIES AGAIN
% [ Vp_est1,Vs_est1,Esurf1,vp_range,vs_range] = ...
%         surf_vp_vs_est([data_stk(:,1,1),data_stk(:,2,1)],tt_stk(:,1),rayp_sdeg2skm(rayp_av(ip)),[-5 5],1 );
% [ Vp_est2,Vs_est2,Esurf2] = ...
%         surf_vp_vs_est([data_stk(:,1,2),data_stk(:,2,2)],tt_stk(:,2),rayp_sdeg2skm(rayp_av(ip)),[-3 4],1 );
% [~,x1,y1] = mingrid(2*Esurf1/mingrid(Esurf1)+Esurf2/mingrid(Esurf2));
% vs_range(x1)
% vp_range(y1)

%% save
if ifPSV
    components = {'P','SV','SH'}; 
    datatype = 'dataPSVSH';
else
    components = eqar.components;  
    datatype = 'dataZRT';
end

if size(plot_dat,2)<2 || isempty(plot_dat{8,2})
    phkeep = phases(1); nevinds = nevinds(1); evinds_use = evinds_use(1); rayp_av = rayp_av(1); gcarc_av = gcarc_av(1); seaz_av = seaz_av(1); data_stk = data_stk(:,:,1); tt_stk = tt_stk(:,1);
elseif size(plot_dat,2)>1 && isempty(plot_dat{8,1})
    phkeep = phases(2); nevinds = nevinds(2); evinds_use = evinds_use(2); rayp_av = rayp_av(2); gcarc_av = gcarc_av(2); seaz_av = seaz_av(2); data_stk = data_stk(:,:,2); tt_stk = tt_stk(:,2);
    plot_dat(:,1) = plot_dat(:,2);
else
    phkeep = phases;
end
plot_dat{1,3} = phkeep;
plot_dat{2,3} = components;

% Figure5

avar = struct('sta',eqar.sta,'nwk',eqar.nwk,'phases',{phkeep},'components',{components},...
              'slat',eqar.slat,'slon',eqar.slon,'selev',eqar.selev,...
              'norids',eqar.norids,'evinds',{evinds_use},'nevinds',nevinds,...
              'elats',eqar.elats,'elons',eqar.elons,'edeps',eqar.edeps,...
              'emags',eqar.emags,'evtimes',{eqar.evtimes},...
              'seaz',seaz_av,'seazs',eqar.seazs,...
              'gcarc',gcarc_av,'gcarcs',eqar.gcarcs,...
              'rayp',rayp_av,'rayps',eqar.rayps,...
              'SNRmin',SNRmin,'acormin',acormin,...
              'xcor_filt',xcor_filtfs,'xcor_win',xcor_win,...
              'filtfs',filtfs,'SNR_stk',SNR_stk,'Vp_Vs_surf',[Vp_surf,Vs_surf],...
              datatype,data_stk,'tt',tt_stk,'stackmethod','doublestack');
          
if ifsave
    if ~any(avar.nevinds>=min4stack), continue; end
    
    fprintf('SAVING----------\n')
    
    arfile = sprintf('avar_dat_%s_%s_%.0f_%.0f',station,network,mean(gcarc_av),mean(seaz_av));
    if ifpatricide, arfile = [arfile,'_patricide']; end 
    save(arfile,'avar');
    
    plotdatfile = sprintf('plotting_dat_%s_%s_%.0f_%.0f',station,network,mean(gcarc_av),mean(seaz_av));
    if ifpatricide, plotdatfile = [plotdatfile,'_patricide']; end
    save(plotdatfile,'plot_dat')
    
    if exist(odir,'dir'); movefile([arfile,'.mat'],odir); end
    if exist(odir,'dir'); movefile([plotdatfile,'.mat'],odir); end
end

fprintf('Cluster %.0f complete\n\n',icl)
end % on cluster loop


end % station loop
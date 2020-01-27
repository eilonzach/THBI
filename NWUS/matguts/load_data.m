function [ trudata,zeroDstr ] = load_data( par )
%[ trudata,cheatstr ] = load_data( par )

if nargin < 6
    baz = [];
end

P_win = [-5 25];
S_win = [-35 5] ;
tapertime = 2;
samprate = 10;
zeroDstr = []; % whether or not the daughter component is zeroed at the time of the parent
% zeroDstr = '_zeroDstr';

wd = pwd;

BWavardir = par.data.avardir;
sta = par.data.stadeets.sta;
nwk = par.data.stadeets.nwk;

%% work out which data types to grab
allpdytp = parse_dtype_all(par);

%% ------------------------------------------------------------------
%% --------------------  BODY WAVE STACKED DATA  --------------------
if any(strcmp(allpdytp(:,1),'BW'))
    cd(BWavardir);
    %% Seismogram data
    datfiles = dir(sprintf('avar_dat_%s_%s_*%s.mat',sta,nwk,zeroDstr));

    if isempty(datfiles)
        trudata=[];
        return
    end

    clear('allavars','seazmean');
    for id = 1:length(datfiles)
        load(datfiles(id).name);
        allavars(id) = avar;
        seazmean(id) = mean(avar.seaz);
    end
    % find clusters of backazimuths
    clusts = eqcluster(seazmean,90*ones(size(seazmean)),0.5,30,10,0);

    np = 0;
    ns = 0;
    for ii = 1:length(allavars)
        avar = allavars(ii);
        clust = clusts(ii);
        Pind = strcmp(avar.phases,'P') | strcmp(avar.phases,'Ps') ;
        Sind = strcmp(avar.phases,'S') | strcmp(avar.phases,'Sp') ;

        %% BW_Ps ==> flip Z to 'up', taper, downsample, window

        if any(Pind)
            Pdat = avar.dataPSVSH(:,:,Pind); 
            Pdat_t  = flat_hanning_win(avar.tt(:,Pind),Pdat,P_win(1),P_win(2),tapertime);
            Pdat_td = downsamp(Pdat_t,unique(round(1./diff(avar.tt(:,Pind)))),samprate);
            tt_d = avar.tt(1,Pind) + [0:size(Pdat_td,1)-1]'/samprate;
            tt_w = P_win(1) + [0:(diff(P_win)*samprate-1)]'./samprate;
            Pdat_tdw = interp1(tt_d,Pdat_td,tt_w);
            np = np+1;
            BW_Ps(np,1) = struct('PSV',Pdat_tdw(:,1:2),'tt',tt_w,...
                                 'clust',clust,'rayp',avar.rayp(Pind),'seaz',avar.seaz(Pind),'gcarc',avar.gcarc(Pind),...
                                 'samprate',samprate,'nsamp',size(Pdat_tdw,1),...
                                 'Vp_surf',avar.Vp_Vs_surf(1),'Vs_surf',avar.Vp_Vs_surf(2));
        end

        %% BW_Sp ==> flip Z to 'up', taper, downsample, window
        if any(Sind)
            Sdat = avar.dataPSVSH(:,:,Sind); 
            Sdat_t  = flat_hanning_win(avar.tt(:,Sind),Sdat,S_win(1),S_win(2),tapertime);
            Sdat_td = downsamp(Sdat_t,unique(round(1./diff(avar.tt(:,Sind)))),samprate);
            tt_d = avar.tt(1,Sind) + [0:size(Sdat_td,1)-1]'/samprate;
            tt_w = S_win(1) + [0:(diff(S_win)*samprate-1)]'./samprate;
            Sdat_tdw = interp1(tt_d,Sdat_td,tt_w);
            ns = ns+1;
            BW_Sp(ns,1) = struct('PSV',Sdat_tdw(:,1:2),'tt',tt_w,...
                                 'clust',clust,'rayp',avar.rayp(Sind),'seaz',avar.seaz(Sind),'gcarc',avar.gcarc(Sind),...
                                 'samprate',samprate,'nsamp',size(Sdat_tdw,1),...
                                 'Vp_surf',avar.Vp_Vs_surf(1),'Vs_surf',avar.Vp_Vs_surf(2));
        end

    % ------------------------- OLD - BASED ON ZRT  ---------------------------
    %         if any(Pind)
    %             Pdat = avar.dataZRT(:,:,Pind); 
    %             Pdat(:,1) = -Pdat(:,1);
    %             Pdat_t  = flat_hanning_win(avar.tt(:,Pind),Pdat,P_win(1),P_win(2),tapertime);
    %             Pdat_td = downsamp(Pdat_t,unique(round(1./diff(avar.tt(:,Pind)))),samprate);
    %             tt_d = avar.tt(1,Pind) + [0:size(Pdat_td,1)-1]'/samprate;
    %             tt_w = P_win(1) + [0:(diff(P_win)*samprate-1)]'./samprate;
    %             Pdat_tdw = interp1(tt_d,Pdat_td,tt_w);
    %             % flip polarity based on Zmax (=> main arrival positive for all)
    %             if maxab(Pdat_tdw(:,1))<0, Pdat_tdw = -Pdat_tdw; end
    %             np = np+1;
    %             BW_Ps(np,1) = struct('ZRT',Pdat_tdw,'tt',tt_w,'rayp',avar.rayp(Pind),'samprate',samprate,'nsamp',size(Pdat_tdw,1));
    %         end
    %         
    %         %% SpRF ==> flip Z to 'up', taper, downsample, window
    %         if any(Sind)
    %             Sdat = avar.dataZRT(:,:,Sind); 
    %             Sdat(:,1) = -Sdat(:,1); 
    %             Sdat_t  = flat_hanning_win(avar.tt(:,Sind),Sdat,S_win(1),S_win(2),tapertime);
    %             Sdat_td = downsamp(Sdat_t,unique(round(1./diff(avar.tt(:,Sind)))),samprate);
    %             tt_d = avar.tt(1,Sind) + [0:size(Sdat_td,1)-1]'/samprate;
    %             tt_w = S_win(1) + [0:(diff(S_win)*samprate-1)]'./samprate;
    %             Sdat_tdw = interp1(tt_d,Sdat_td,tt_w);
    %             % flip polarity based on Rmax (=> main arrival positive for all)
    %             if maxab(Sdat_tdw(:,2))<0, Sdat_tdw = -Sdat_tdw; end
    %             ns = ns+1;
    %             BW_Sp(ns,1) = struct('ZRT',Sdat_tdw,'tt',tt_w,'rayp',avar.rayp(Sind),'samprate',samprate,'nsamp',size(Sdat_tdw,1));
    %         end
    % ------------------------- OLD - BASED ON ZRT  ---------------------------
    end % end loop on data files

    cd(wd);
end

%% ------------------------------------------------------------------
%% --------------------  RECEIVER FUNCTION DATA  --------------------
if any(strcmp(allpdytp(:,1),'RF'))
    for idt = find(strcmp(allpdytp(:,1),'RF'))
    % should we be pulling in CCP data?
        if strcmpi(allpdytp(idt,3),'CCP')
            fprintf('\n Extracting %s "RF" from CCP stack\n',allpdytp{idt,2})
            % addpath to CCP data
            addpath('/Volumes/data/models_seismic/');
            slat = par.data.stadeets.Latitude;
            slon = par.data.stadeets.Longitude;

            % grab RF(z) from CCP
            % negative P polarity means positive velocity gradient (e.g.
            % moho), consistent with the forward models
            [RF_dz,zz_d] = load_CCP_RFz( [slat,slon],0 ); 
            zz_d = double(zz_d);
            dz = unique(diff(zz_d));
            
            % window to desired depth range
            Zwin = par.datprocess.CCP.Zwin.def;
            taperz = par.datprocess.CCP.taperz;
            RF_dzw = flat_hanning_win(zz_d,RF_dz,Zwin(1)-taperz/2,Zwin(2)+taperz/2,taperz);
            
            

            % migrate to time
%             if strcmpi(par.datprocess.CCP.migratemodel,'PREM')
%                 premmod = prem;
%                 model = struct('z',linterp(premmod.depth(1:20),premmod.depth(1:20),zz(zz<280)),...
%                                'VS',linterp(premmod.depth(1:20),premmod.vs(1:20),zz(zz<280)),...
%                                'VP',linterp(premmod.depth(1:20),premmod.vp(1:20),zz(zz<280)));
%                 RFzz = interp1(zz(zz<280),RFz(zz<280),model.z);
%             end        
%             rayp = rayp_sdeg2skm(par.datprocess.CCP.rayp_S); % N.B. should this be at a depth of zero, or the bottom of the model?
%             fprintf(' Migrating RF(z) to RF(t) using %s, rayp = %.3f s/deg\n',par.datprocess.CCP.migratemodel,par.datprocess.CCP.rayp_S)
%             [RF_tmig,tt_tmig] = migrate_depth2time(RFzz,model,rayp,allpdytp{idt,2},1./samprate);
%             % convert migrated RF to appropriate length defined by Twin
%             tt_RF = [Twin(1):1./samprate:Twin(2) - 1./samprate]';
%             RF_d  = interp1(tt_tmig,RF_tmig,tt_RF,'linear',0);


            % make synthetic parent
            RF_p = synthtrace(par.datprocess.CCP.parent_zw/2 + max(Zwin),...
                              par.datprocess.CCP.parent_zw,...
                              1,dz,'gauss',0);
            zz = [-par.datprocess.CCP.parent_zw/2:dz:max(Zwin)-dz]';
            RF_d = interp1(zz_d,RF_dzw,zz,'linear',0);
                                      % 

%             fprintf('NOTE WHAT TO DO ABOUT SURFV FOR CCP??\n');
            RF_struct = struct('PSV',[RF_d,RF_p],'zz',zz,...
                 'rayp',par.datprocess.CCP.rayp_S,...
                 'dz',dz,'nsamp',size(zz,1),...
                 'dof_per_z',1/par.datprocess.CCP.parent_zw); % say parent pulse width is 1 dof
            eval(sprintf('RF_%s_ccp = RF_struct;',allpdytp{idt,2}));
        else
            error('expecting CCP...')   ; 
        end
    end
end

%% ------------------------------------------------------------------
%% ----------------------  SURFACE WAVE DATA  ----------------------
if any(strcmp(allpdytp(:,1),'SW'))

    %% Phase velocity data
    seismoddir = '/Volumes/data/models_seismic/';
    if ~exist(seismoddir,'dir'), seismoddir = regexprep(seismoddir,'~','/Volumes/zeilon'); end 
    addpath(seismoddir);

    slat = par.data.stadeets.Latitude;
    slon = par.data.stadeets.Longitude;

    %% -------- Rayleigh waves
    if any(strcmp(allpdytp(strcmp(allpdytp(:,1),'SW'),2),'Ray'))
        [Rperiods,RphV]  = Rph_dispcurve_latlon( slat,slon); % grab composite of AN and EQ
        % err = error_curve_EQ_latlon( 1./periods,avar.slat,avar.slon,phVerrordir);

        if ~isempty(RphV)
            [Rperiods,iT] = sort(Rperiods);
            RphV = RphV(iT);
            SW_Ray_phV = struct('periods',Rperiods,'phV',RphV,'sigma',[]);
        else
            SW_Ray_phV=[];
        end
    end

    %% -------- Love waves
	if any(strcmp(allpdytp(strcmp(allpdytp(:,1),'SW'),2),'Lov'))
        [Lperiods,LphV]  = Lph_dispcurve_latlon( slat,slon); % grab composite of AN and EQ

        if ~isempty(LphV)
            [Lperiods,iT] = sort(Lperiods);
            LphV = LphV(iT);
            SW_Lov_phV = struct('periods',Lperiods,'phV',LphV,'sigma',[]);
        else
            SW_Lov_phV=[];
        end
    end

    %% -------- Rayleigh HV ratios
	if any(strcmp(allpdytp(strcmp(allpdytp(:,1),'SW'),2),'HV'))
        [ HVratios,HVstds,HVperiods] = disp_curve_HV( round_level([slat,slon],0.1),0,[],10,1.5);

        [HVperiods,iT] = sort(HVperiods);
        HVratios = HVratios(iT);
        HVstds = HVstds(iT);

        gdT = ~isnan(HVratios);
        HVperiods = HVperiods(gdT);
        HVratios = HVratios(gdT);
        HVstds = HVstds(gdT);

        if ~isempty(HVratios)
            SW_HV = struct('periods',HVperiods,'HVr',HVratios,'sigma',HVstds);
        else
            SW_HV=[];
        end
    end

end

%% ------------------------------------------------------------------
%% --------------------  H-K stack  --------------------
if any(strcmp(allpdytp(:,1),'HKstack'))
	fprintf('\n pulling out EARS H-K stack\n')
    EARSdir = '/Volumes/data/models_seismic/US_EARS';
    load(sprintf('%s/EARS_HKStack_%s_%s.mat',EARSdir,nwk,sta));
    if isempty(hkstack.Nobs)
        fprintf('\t Not sure how many EQ in EARS obs - assigning 100\n');
        hkstack.Nobs = 100;
    end
    HKstack_P = hkstack;
    HKstack_P.Esum = HKstack_P.Esum-mingrid(HKstack_P.Esum);
end



%% collate
if ~exist('BW_Ps','var'),       BW_Ps = struct([]); end
if ~exist('BW_Sp','var'),       BW_Sp = struct([]); end
if ~exist('SW_Ray_phV','var'),  SW_Ray_phV = []; end
if ~exist('SW_Lov_phV','var'),  SW_Lov_phV = []; end
if ~exist('SW_HV','var'),       SW_HV = []; end
if ~exist('RF_Ps','var'),       RF_Ps = struct([]); end
if ~exist('RF_Sp','var'),       RF_Sp = struct([]); end
if ~exist('RF_Ps_ccp','var'),   RF_Ps_ccp = struct([]); end
if ~exist('RF_Sp_ccp','var'),   RF_Sp_ccp = struct([]); end
if ~exist('HKstack_P','var'),   HKstack_P = struct([]); end


trudata = struct('BW_Ps',BW_Ps,'BW_Sp',BW_Sp,...
                 'RF_Ps',RF_Ps,'RF_Sp',RF_Sp,...
                 'RF_Ps_ccp',RF_Ps_ccp,'RF_Sp_ccp',RF_Sp_ccp,...
                 'HKstack_P',HKstack_P,...
                 'SW_Ray_phV',SW_Ray_phV,'SW_Lov_phV',SW_Lov_phV,'SW_HV',SW_HV);


cd(wd);


function [ trudata,cheatstr ] = load_data( avardir,sta,nwk,gc )
%[ trudata,cheatstr ] = load_data( avardir,sta,nwk,gc )

if nargin < 6
    baz = [];
end

P_win = [-5 35];
S_win = [-35 5] ;
tapertime = 2;
samprate = 10;
cheatstr = [];
% cheatstr = '_cheat';

wd = pwd;

%     addpath('matguts')
cd(avardir);
%% Sesismogram data
datfiles = {};
for gcmid = gc(:)'
for gcwind = -4:4;
    df = dir(sprintf('avar_dat_%s_%s_%.0f_*%s.mat',sta,nwk,gcmid+gcwind,cheatstr));
    datfiles = [datfiles;{df.name}']; %#ok<AGROW>
end
end
datfiles = unique(datfiles);


np = 0;
ns = 0;
for ii = 1:length(datfiles)
    load(datfiles{ii});
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
        BW_Ps(np,1) = struct('PSV',Pdat_tdw(:,1:2),'tt',tt_w,'rayp',avar.rayp(Pind),'samprate',samprate,'nsamp',size(Pdat_tdw,1),'Vp_surf',avar.Vp_Vs_surf(1),'Vs_surf',avar.Vp_Vs_surf(2));
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
        BW_Sp(ns,1) = struct('PSV',Sdat_tdw(:,1:2),'tt',tt_w,'rayp',avar.rayp(Sind),'samprate',samprate,'nsamp',size(Sdat_tdw,1),'Vp_surf',avar.Vp_Vs_surf(1),'Vs_surf',avar.Vp_Vs_surf(2));
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

%% Phase velocity data
seismoddir = '~/Work/data/models_seismic/';
addpath(seismoddir);

% -------- Rayleigh waves
[Rperiods,RphV]  = Rph_dispcurve_latlon( avar.slat,avar.slon); % grab composite of AN and EQ
% err = error_curve_EQ_latlon( 1./periods,avar.slat,avar.slon,phVerrordir);

[Rperiods,iT] = sort(Rperiods);
RphV = RphV(iT);

SW_Ray = struct('periods',Rperiods,'phV',RphV,'sigma',0.01);

% -------- Love waves
[Lperiods,LphV]  = Lph_dispcurve_latlon( avar.slat,avar.slon); % grab composite of AN and EQ

[Lperiods,iT] = sort(Lperiods);
LphV = LphV(iT);

SW_Lov = struct('periods',Lperiods,'phV',LphV,'sigma',0.01);



%% collate
trudata = struct('BW_Ps',BW_Ps,'BW_Sp',BW_Sp,'SW_Ray',SW_Ray,'SW_Lov',SW_Lov);


cd(wd);

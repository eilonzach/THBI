function [ trudata ] = load_data( project,sta,nwk,gc )
%[ trudata ] = load_data( project,sta,nwk,gc )

if nargin < 6
    baz = [];
end

P_win = [-5 35];
S_win = [-35 5] ;
tapertime = 2;
samprate = 10;

if strcmp(project,'WYOMING');
    wd = '~/Documents/MATLAB/BayesianJointInv/WYOMING/';
    addpath('matguts')
    cd(wd);
    
    cd('AVARS');
    %% Sesismogram data
    datfiles = {};
    for gcwind = -4:4;
    	df = dir(sprintf('avar_dat_%s_%s_%.0f_*.mat',sta,nwk,gc+gcwind));
        datfiles = [datfiles;{df.name}']; %#ok<AGROW>
    end
    
    np = 0;
    ns = 0;
    for ii = 1:length(datfiles)
        load(datfiles{ii});
        Pind = strcmp(avar.phases,'P');
        Sind = strcmp(avar.phases,'S');
        
        %% PsRF ==> flip Z to 'up', taper, downsample, window
        if any(Pind)
            Pdat = avar.dataZRT(:,:,Pind); 
            Pdat(:,1) = -Pdat(:,1);
            Pdat_t  = flat_hanning_win(avar.tt(:,Pind),Pdat,P_win(1),P_win(2),tapertime);
            Pdat_td = downsamp(Pdat_t,unique(round(1./diff(avar.tt(:,Pind)))),samprate);
            tt_d = avar.tt(1,Pind) + [0:size(Pdat_td,1)-1]'/samprate;
            tt_w = P_win(1) + [0:(diff(P_win)*samprate-1)]'./samprate;
            Pdat_tdw = interp1(tt_d,Pdat_td,tt_w);
            % flip polarity based on Zmax (=> main arrival positive for all)
            if maxab(Pdat_tdw(:,1))<0, Pdat_tdw = -Pdat_tdw; end
            np = np+1;
            PsRF(np,1) = struct('ZRT',Pdat_tdw,'tt',tt_w,'rayp',avar.rayp(Pind),'samprate',samprate,'nsamp',size(Pdat_tdw,1));
        end
        
        %% SpRF ==> flip Z to 'up', taper, downsample, window
        if any(Sind)
            Sdat = avar.dataZRT(:,:,Sind); 
            Sdat(:,1) = -Sdat(:,1); % NOW WE DO FLIP Z... previously WHY NOT FLIP Z???
            Sdat_t  = flat_hanning_win(avar.tt(:,Sind),Sdat,S_win(1),S_win(2),tapertime);
            Sdat_td = downsamp(Sdat_t,unique(round(1./diff(avar.tt(:,Sind)))),samprate);
            tt_d = avar.tt(1,Sind) + [0:size(Sdat_td,1)-1]'/samprate;
            tt_w = S_win(1) + [0:(diff(S_win)*samprate-1)]'./samprate;
            Sdat_tdw = interp1(tt_d,Sdat_td,tt_w);
            % flip polarity based on Rmax (=> main arrival positive for all)
            if maxab(Sdat_tdw(:,2))<0, Sdat_tdw = -Sdat_tdw; end
            ns = ns+1;
            SpRF(ns,1) = struct('ZRT',Sdat_tdw,'tt',tt_w,'rayp',avar.rayp(Sind),'samprate',samprate,'nsamp',size(Sdat_tdw,1));
        end
    end % end loop on data files
        
    cd(wd);
    
    %% Phase velocity data
    phVdir = '~/Dropbox/Dave_Li_phV/2D_phase_velocities/'; % need final slash
    phVerrordir = '~/Dropbox/Dave_Li_phV/Errors/'; % need final slash

    [ freqs ] = get_freqs(phVdir);

    [periods,phV]  = disp_curve_latlon( avar.slat,avar.slon); % grab composite of AN and EQ
    err = error_curve_EQ_latlon( freqs,avar.slat,avar.slon,phVerrordir);
    
    [periods,iT] = sort(periods);
    phV = phV(iT);
    
    SW = struct('periods',periods,'phV',phV,'sigma',mean(err));
    
    %% collate
    trudata = struct('PsRF',PsRF,'SpRF',SpRF,'SW',SW);

end


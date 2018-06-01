% For all of the evids used in each avar stack, go back to the original
% data, calculate individual receiver functions (migrated to common Vs
% profile) and stack these. Compare this stack of RFs to both the RF of the
% stacked data and to the RF of the final modelled data. 
% clear

avardir = 'STA_inversions/ECSD_dat20/';% need final slash
datadir = 'DATA/';% need final slash
station = 'ECSD';
network = 'US';
gcarc_av = 67;
baz_av = 155;

RFphase = {'Ps','Sp'};
dat_win = [-40 40]; % should have main phase roughly centred
dat_filtfs = [[0.02;100],[0.02;100]]; % p,s [flo fhi]
npoles = 4;
dat_win = [-30 30];

% Receiver function
TB = 1.5;%1.5
NT = 2;%2
w_len = diff(dat_win);
Poverlap = 0.9;


ifsave = true;

%% load things
% raw data
% load(sprintf('%sdat_%s_%s_30to75.mat',datadir,station,network));

% avar
load(sprintf('%savar_dat_%s_%s_%.0f_%.0f_cheat.mat',avardir,station,network,gcarc_av,baz_av));

addpath('plotting','matguts','../functions')
addpath('~/Dropbox/MATLAB/myFUNCTIONS/','~/Dropbox/MATLAB/seis_tools/')

for ip = 1:length(RFphase)
    ph = RFphase{ip};
    ipha = find(strcmp(avar.phases,ph)); if isempty(ipha), continue, end
    iphd = find(strcmp(eqar.phases,ph(1)));
    RFs = nan(avar.nevinds(ipha),1);
    bazs = nan(avar.nevinds(ipha),1);

%% make individual RFs
for ie = 1:avar.nevinds(ipha)
    try evind = avar.evinds{ipha}(ie); catch, evind = avar.evinds(ie); end
    rayp_skm = rayp_sdeg2skm(eqar.rayps(evind,iphd));
    bazs(ie) = avar.seazs(evind);
    
    %% Get data, filter
    datZ = eqar.dataZRT(:,evind,1,iphd);
    datR = eqar.dataZRT(:,evind,2,iphd);
    tt = eqar.tt(:,evind,iphd); samprate = round(1./diff(tt(1:2)));

    % rotate to P-SV
    [datP,datSV] = Rotate_XZ_to_PSV(datR,datZ,avar.Vp_Vs_surf(1),avar.Vp_Vs_surf(2),rayp_skm);

    % filter, clean, window
    cp = struct('samprate',samprate,'fhi',dat_filtfs(2,ip),'flo',dat_filtfs(1,ip),...
            'pretime',-tt(1),'prex',-dat_win(1),'postx',dat_win(2),...
            'taperx',1./diff(dat_win),'npoles',npoles,'norm',1);
    [ datP,~,~,~,~,~] = data_clean(  datP,cp ); 
    [ datSV,~,~,~,~,tt] = data_clean(  datSV,cp ); 

    if strcmp(ph,'Sp')
        Parent = datSV; Daughter = datP;
    elseif strcmp(ph,'Ps')
        Parent = datP; Daughter = datSV;
    end

    % taper Parent and daughter
    Parent = flat_hanning_win(tt,Parent,-5,5,1);
    if strcmp(ph,'Sp')
        Daughter = flat_hanning_win(tt,Daughter,dat_win(1),5,1);
    elseif strcmp(ph,'Ps')
        Daughter = flat_hanning_win(tt,Daughter,-5,dat_win(2),1);
    end

    % sign of time direction
    if strcmp(ph,'Ps'), phsgn = 1; elseif strcmp(ph,'Sp'), phsgn = -1; end


    [tt_RF, RF] = ETMTM(Parent',Daughter',ph,TB,NT,'synth',1./samprate,w_len,Poverlap);
    if size(RFs,2)==1; RFs = nan(eqar.norids,length(tt_RF)); end
    RFs(ie,:) = RF(:)'./max(abs(RF));

end

%% Make stack of RFs
RFsum = nansum(RFs)/sum(all(~isnan(RFs')));

%% Make RF of stack
    % filter, clean, window
    cp = struct('samprate',samprate,'fhi',dat_filtfs(2,ip),'flo',dat_filtfs(1,ip),...
            'pretime',-avar.tt(1,ipha),'prex',-dat_win(1),'postx',dat_win(2),...
            'taperx',1./diff(dat_win),'npoles',npoles,'norm',1);
    [ datP,~,~,~,~,~]   = data_clean( avar.dataPSVSH(:,1,ipha),cp ); 
    [ datSV,~,~,~,~,tt] = data_clean( avar.dataPSVSH(:,2,ipha),cp ); 

    if strcmp(ph,'Sp')
        Parent = datSV; Daughter = datP;
    elseif strcmp(ph,'Ps')
        Parent = datP; Daughter = datSV;
    end

    % taper Parent and daughter
    Parent = flat_hanning_win(tt,Parent,-5,5,1);
    if strcmp(ph,'Sp')
        Daughter = flat_hanning_win(tt,Daughter,dat_win(1),5,1);
    elseif strcmp(ph,'Ps')
        Daughter = flat_hanning_win(tt,Daughter,-5,dat_win(2),1);
    end

    % sign of time direction
    if strcmp(ph,'Ps'), phsgn = 1; elseif strcmp(ph,'Sp'), phsgn = -1; end

    [tt_RFstk, RFstk] = ETMTM(Parent',Daughter',ph,TB,NT,'synth',1./samprate,w_len,Poverlap);







%% PLOT
figure(86+ip),clf, set(gcf,'pos',[384 139 1102 860]);hold on
% RFs_filn = RFs; RFs_filn(RFs_filn<0.2)=nan;
for ie = 1:avar.nevinds(ipha)
    plot_RF( gca,0.5*RFs(ie,:),phsgn*tt_RF,1*rms(RFs(ie,:)),bazs(ie) );
end
dx = max(bazs)-avar.seaz(ipha); RFsx = avar.seaz(ipha)+dx+4; sRFx = avar.seaz(ipha)+dx+2; 

plot_RF( gca,2*RFsum,phsgn*tt_RF,2*rms(RFsum),sRFx );
plot_RF( gca,2*RFstk,phsgn*tt_RFstk,2*rms(RFstk),RFsx );
plot(sRFx + 2*RFsum,phsgn*tt_RF,'k','linewidth',1.5)
plot(RFsx + 2*RFstk,phsgn*tt_RFstk,'k','linewidth',1.5)

set(gca,'ylim',[-2 30],'ydir','reverse','xlim',[min(bazs)-1 max(bazs)+6],...
    'xtick',[floor(min(bazs)):1:ceil(max(bazs))]);


text(sRFx+0.5,22,'Stack of RFs','fontweight','bold','rotation',270,'fontsize',18)
text(RFsx+0.5,22,'RF of stack','fontweight','bold','rotation',270,'fontsize',18)

title(sprintf('%s  %s  receiver functions ($\\Delta$%.0f $\\phi$%.0f)',station,ph,gcarc_av,baz_av),'fontsize',24,'interpreter','latex')

save2jpg(86+ip,sprintf('%sRF_%.0f_%.0f',ph,gcarc_av,baz_av),avardir);

end








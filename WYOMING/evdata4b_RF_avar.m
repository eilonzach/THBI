close all
clear all
avardir = 'STA_inversions/RSSD_dat20/';% need final slash
station = 'RSSD';
network = 'IU';
gcarc_av = 69;
baz_av = 147;

RFphase = {'Ps','Sp'};
filtfs = [[1/50 ;1/4],[1/2000;1/5]];
npoles = 4;
wdo = [-40 40];
tapertime = 5;

ifsave = false;

load(sprintf('%savar_dat_%s_%s_%.0f_%.0f',avardir,station,network,gcarc_av,baz_av));

addpath('plotting','matguts','../functions')
% [~,~,vslay] = RD_1D_Vprofile;
% Vs_surf = vslay(1);
% Vp_surf = sed_vs2vp(Vs_surf);
% Vs_surf = 3.5;
% Vp_surf = 1.68*Vs_surf;



for ip = 1:length(RFphase)
    ph = RFphase{ip};
    
    rayp_skm = rayp_sdeg2skm(avar.rayp(ip));
    
    %% Get data, filter
    datSV = avar.dataPSVSH(:,2,ip);
    datP = avar.dataPSVSH(:,1,ip);
    tt = avar.tt(:,ip); 
    
    %% Window
    datSV = flat_hanning_win(tt,datSV,wdo(1),wdo(2),tapertime);
    datP = flat_hanning_win(tt,datP,wdo(1),wdo(2),tapertime);
    datSV = datSV(tt>=wdo(1) & tt < wdo(2));
    datP = datP(tt>=wdo(1) & tt < wdo(2));
    tt = tt(tt>=wdo(1) & tt < wdo(2));
    
    figure(ip); 
    subplot(311)
    plot(tt,datSV,'r',tt,datP,'k')
    samprate = unique(round(1./diff(tt)));
    
    datSV = filt_quick(datSV,filtfs(1,ip),filtfs(2,ip),1./samprate,npoles);
    datP = filt_quick(datP,filtfs(1,ip),filtfs(2,ip),1./samprate,npoles);
    
    subplot(312)
    plot(tt,datSV,'r',tt,datP,'k'),%pause
    

    % rotate to P-SV
%     [P,SV] = Rotate_XZ_to_PSV(datSV,datP,Vp_surf,Vs_surf,rayp_skm);
    
    subplot(313)
    plot(tt,datP,'r',tt,datSV,'k'),%pause
    
    pow = max(max(abs([datP,datSV])));
    datP = datP./pow; 
    datSV = datSV./pow;
    plot(tt,datP,'k',tt,datSV,'r')

    % Receiver function
    TB = 1.5;%1.5
    NT = 2;%2
    w_len = 60;
    Poverlap = 0.9;
    if strcmp(ph,'Sp')
        Parent = datSV; Daughter = datP;
    elseif strcmp(ph,'Ps')
        Parent = datP; Daughter = datSV;
    end

    % taper Parent and daughter
    % Parent = flat_hanning_win(tt,Parent,-10,10,1);
    if strcmp(ph,'Sp')
        Daughter = flat_hanning_win(tt,Daughter,wdo(1),5,1);
        Parent   = flat_hanning_win(tt,Parent,wdo(1),5,1);
    elseif strcmp(ph,'Ps')
        Daughter = flat_hanning_win(tt,Daughter,-5,wdo(2),1);
        Parent   = flat_hanning_win(tt,Parent,-5,wdo(2),1);
    end

    % sign of time direction
    if strcmp(ph,'Ps'), phsgn = 1; elseif strcmp(ph,'Sp'), phsgn = -1; end


    [tt_RF, RF] = ETMTM(Parent',Daughter',ph,TB,NT,'synth',1./samprate,w_len,Poverlap);
    RF = RF(:)'./max(abs(RF));
    % % RF
    % figure(68), clf; set(gcf,'pos',[403 176 1315 529])
    % ax1 = subplot(311); hold on, set(gca,'fontsize',15)
    % plot(tt,P,'linewidth',2); 
    % title(ph,'fontsize',25), 
    %     
    % ax2=subplot(312); hold on, set(gca,'fontsize',15)
    % plot(tt,SV,'linewidth',2), 
    % 
    % ax3=subplot(313); hold on, set(gca,'fontsize',15)
    % plot(ax3,tt_RF,RF,'linewidth',2), 
    % 
    % ylabs = {'P','SV','RF'};
    % ii=0;
    % for axii = [ax1,ax2,ax3]
    %     ii = ii+1;
    %     set(axii,'xlim',sort(phsgn*[-2 30])); 
    %     if strcmp(ph,'Sp'), set(axii,'xdir','reverse'); end
    %     hy3 = ylabel(axii,['\textbf{',ylabs{ii},'}'],'fontsize',22,'interpreter','latex','rotation',0,'verticalalignment','middle');
    %     set(hy3,'position',[-3.5*phsgn,0,1]);
    % end
    % set(ax3,'ylim',max(abs(RF))*1.2*[-1 1])
    % % % time label
    % hx = xlabel(ax3,'\textbf{Time (s)}  $\mathbf{\longrightarrow}$','fontsize',22,'interpreter','latex');
    % if strcmp(ph,'Sp'), 
    %     set(hx,'string','$\mathbf{\longleftarrow}$  \textbf{Time (s)}','fontweight','bold')
    % end
    % 
    % pause

    figure(87),
    if ip==1, clf, set(gcf,'pos',[384 139 600 860]); end

    ax = subplot(1,2,ip); hold on

    plot_RF( ax,phsgn*10*RF,phsgn*tt_RF,20*rms(RF));
    plot(ax,phsgn*10*RF,phsgn*tt_RF,'k','linewidth',1.5)

    set(ax,'ylim',sort([-2 30]),'ydir','reverse','fontsize',15);
    title(ax,RFphase{ip},'fontsize',20)

    if ip==1, ylabel('Time (seconds after/before parent)','fontsize',20); end
    if ip==1, xlabel('RF amplitude','fontsize',20); end
    if ip==2, xlabel('-RF amplitude','fontsize',20); end



end

if ifsave 
    save2pdf(87,sprintf('RF_avar_%s_%s_%.0f',station, network, gcarc_av),'figs')
end







clear all
station = 'LKWY';
network = 'US';

RFphase = {'Ps','Sp'};
filtfs = [[0.001;0.5],[0.001;0.5]];
npoles = 4;
wdo = [-30 30];

SNRmin = 3;

ifsave = false;

load(sprintf('dat_%s_%s',station,network));

addpath('plotting','matguts','../functions')
[~,~,vslay] = RD_1D_Vprofile;
Vs_surf = vslay(1);
Vp_surf = sed_vs2vp(Vs_surf);

for ip = 2:length(RFphase)
    ph = RFphase{ip};
    RFs = nan(eqar.norids,1);
for ie = 1:eqar.norids
    if strcmp(ph,'Sp')
        if eqar.SNRs(ie,2) < SNRmin, continue, end
    elseif strcmp(ph,'Ps')
        if eqar.SNRs(ie,1) < SNRmin, continue, end
    end
    
    ie
    rayp_skm = rayp_sdeg2skm(eqar.rayps(ie));
    
    %% Get data, filter
    datR = eqar.dataZRT(:,ie,2,ip);
    datZ = eqar.dataZRT(:,ie,1,ip);
    tt = eqar.tt(:,ie,ip); 
    figure(1); 
    subplot(211)
    plot(tt,datR,'r',tt,datZ,'k')
    samprate = unique(round(1./diff(tt)));
    datR = filt_quick(datR,filtfs(1,ip),filtfs(2,ip),1./samprate,npoles);
    datZ = filt_quick(datZ,filtfs(1,ip),filtfs(2,ip),1./samprate,npoles);
    subplot(212)
    plot(tt,datR,'r',tt,datZ,'k'),%pause

% rotate to P-SV
[P,SV] = Rotate_XZ_to_PSV(datR,datZ,Vp_surf,Vs_surf,rayp_skm);

pow = max(max(abs([P,SV])));
P = P./pow; 
SV = SV./pow;
plot(tt,P,'k',tt,SV,'r')

    
% Receiver function
TB = 1.5;%1.5
NT = 2;%2
w_len = diff(wdo);
Poverlap = 0.9;
if strcmp(ph,'Sp')
    Parent = SV; Daughter = P;
elseif strcmp(ph,'Ps')
    Parent = P; Daughter = SV;
end

% taper Parent and daughter
Parent = flat_hanning_win(tt,Parent,-5,5,1);
if strcmp(ph,'Sp')
    Daughter = flat_hanning_win(tt,Daughter,wdo(1),5,1);
elseif strcmp(ph,'Ps')
    Daughter = flat_hanning_win(tt,Daughter,-5,wdo(2),1);
end

% sign of time direction
if strcmp(ph,'Ps'), phsgn = 1; elseif strcmp(ph,'Sp'), phsgn = -1; end


[tt_RF, RF] = ETMTM(Parent',Daughter',ph,TB,NT,'synth',1./samprate,w_len,Poverlap);
if size(RFs,2)==1; RFs = nan(eqar.norids,length(tt_RF)); end
RFs(ie,:) = RF(:)'./max(abs(RF));
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

end
figure(87),clf, set(gcf,'pos',[384 139 1102 860]);hold on
% RFs_filn = RFs; RFs_filn(RFs_filn<0.2)=nan;
for ie = 1:eqar.norids
%   [h,hline,hpos,hneg] = plot_RF( gca,10*RFs(ie,:),tt_RF,2,eqar.seazs(ie) );
    plot_RF( gca,10*RFs(ie,:),phsgn*tt_RF,20*rms(RFs(ie,:)),eqar.seazs(ie) );
end
% plot(eqar.seazs*ones(1,length(tt_RF)) + 10*RFs,tt_RF,'k','linewidth',0.5)
RFsum = nansum(RFs)/sum(all(~isnan(RFs')));
[h,hline,hpos,hneg] = plot_RF( gca,10*RFsum,phsgn*tt_RF,20*rms(RFsum),400 );
plot(400 + 10*nansum(RFs)/sum(all(~isnan(RFs'))),phsgn*tt_RF,'k','linewidth',1.5)

set(gca,'ylim',sort([-2 30])); 
set(gca,'ydir','reverse');



return
end








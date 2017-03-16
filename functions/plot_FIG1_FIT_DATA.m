function plot_FIG1_FIT_DATA( trudata,savedata,par,ifsave,ofile)
% plot_FIG1_FIT_DATA( trudata,savedata,par,ifsave,ofile)
%   

savedata0 = savedata;

if nargin < 4 || isempty(ifsave)
    ifsave=false;
end

if nargin < 5 || isempty(ofile)
    ofile='fig1_FIT_DATA';
end

ps_xlims = [-1 25];
sp_xlims = [-30 1];
synthcol = 0.7*[1 1 1];

iflodata = false;

downsampx = 5; % 5x down-sample both in time and in every nth trace to plot

figure(61),clf,
if iflodata
    set(gcf,'pos',[015 576 1300 529])
    ax1 = axes('position',[0.05 0.56 0.29 0.40]); hold on
    ax2 = axes('position',[0.05 0.10 0.29 0.40]); hold on
    ax4 = axes('position',[0.36 0.56 0.29 0.40]); hold on
    ax5 = axes('position',[0.36 0.10 0.29 0.40]); hold on
    ax3 = axes('position',[0.71 0.10 0.26 0.84]); hold on
    axs=[ax1,ax2,ax3,ax4,ax5];
else
    set(gcf,'pos',[015 576 1000 529])
    ax1 = axes('position',[0.08 0.56 0.41 0.40]); hold on
    ax2 = axes('position',[0.08 0.10 0.41 0.40]); hold on
    ax3 = axes('position',[0.57 0.10 0.39 0.86]); hold on
    axs=[ax1,ax2,ax3];
end

%%  =========================  PROCESS SAVEDATA  =========================  

if iflodata
    dtypes = {'PsRF','SpRF','PsRF_lo','SpRF_lo'};
    axord = [4,5,1,2,3]
else
    dtypes = {'PsRF','SpRF'};
    axord = [1,2];
end

fprintf('Processing saved raw data... ')
for id = 1:length(dtypes)
    dtype = dtypes{id}; fprintf('%s... ',dtype);
    odat = savedata0.(dtype(1:4));
    odat.tt = round_level(odat.tt(savedata.gdmods(1),:),0.001); if odat.tt(end)==0; odat.tt(end)=nan; end

	%% apply windowing - next will filter, taper
    gdtt = (odat.tt>= par.datprocess.Twin.(dtype)(1)) & (odat.tt < par.datprocess.Twin.(dtype)(2));
    odat.tt = odat.tt(:,gdtt)'; % flip
    odat.Z = odat.Z(savedata.gdmods,gdtt)'; % flip
    odat.R = odat.R(savedata.gdmods,gdtt)'; % flip
    odat.T = odat.T(savedata.gdmods,gdtt)'; % flip

    cp = struct('samprate',trudata.(dtype)(1).samprate,        ...
                'pretime',-par.datprocess.Twin.(dtype)(1),  ...
                'prex',-par.datprocess.Twin.(dtype)(1),     ...
                'postx',par.datprocess.Twin.(dtype)(2),     ...                 
                'fhi',par.datprocess.filtf.(dtype)(1),      ...
                'flo',par.datprocess.filtf.(dtype)(2),      ...
                'taperx',0.06,'npoles',2,'norm',0           );

    %% clean, filter, taper
    odat.Z = data_clean(odat.Z,cp);
    odat.R = data_clean(odat.R,cp);
    odat.T = data_clean(odat.T,cp);

%% normalise to unit energy
if par.datprocess.normdata 
    normf = diag(odat.Z'*odat.Z) + diag(odat.R'*odat.R) +diag(odat.T'*odat.T);
    odat.Z = odat.Z*diag(1./sqrt(normf));
    odat.R = odat.R*diag(1./sqrt(normf));
    odat.T = odat.T*diag(1./sqrt(normf));
end

 
%% decimate
% if par.datprocess.decdata
%     resamprate = cp.fhi*4;
%     tt_new = [odat.tt(1):1./resamprate:dat_out.tt(end)]';
%     dat_out.ZRT = interp1(dat_out.tt,dat_out.ZRT,tt_new);
%     dat_out.tt = tt_new;
%     dat_out.nsamp = length(dat_out.tt);
%     dat_out.samprate=resamprate;
% end

%% spit back out
    savedata.(dtype) = odat;
end
fprintf('\nPlotting... ');

%%  ====================  PLOTTING  ====================

%% RFs
for id = 1:length(dtypes)
    dtype = dtypes{id}; fprintf('%s ',dtype);
    if ~iflodata, if ~isempty(regexp(dtype,'_lo','once')), continue; end, end
    xa = axs(axord(id)); 
    cla(xa); clear('cc_v1','cc_v2')
    samprate = trudata.(dtype)(1).samprate;
    for itr = 1:length(trudata.(dtype))
        fprintf('%.0f ',itr);
        for is = 1:size(savedata.(dtype).Z,2)
            cc_v1(:,is) = conv(trudata.(dtype)(itr).ZRT(:,1),savedata.(dtype).R(:,is),'full'); % Vobs*Hpre
            cc_v2(:,is) = conv(trudata.(dtype)(itr).ZRT(:,2),savedata.(dtype).Z(:,is),'full'); % Hobs*Vpre
        end
        cc_max = max(max(abs([cc_v1;cc_v2])));
        cc_t = [0:size(cc_v1,1)-1]'./samprate;
        
        idspx = [1:downsampx:size(cc_v1,1)]; %DOWNSAMPLE FOR THE PLOT
        idspy = [1:downsampx:size(cc_v1,2)]; %DOWNSAMPLE FOR THE PLOT
        
        plot(xa,cc_t(idspx),cc_v1(idspx,idspy),'linewidth',2.5,'color',synthcol+[0   0.1*(itr-1) 0.2])
        plot(xa,cc_t(idspx),cc_v2(idspx,idspy),'linewidth',1.5,'color',synthcol+[0.2 0.1*(itr-1) 0])
    end
    if strcmp(dtype(1),'P')
        xlims = [0 1.5*trudata.(dtype)(1).nsamp./samprate];
    elseif strcmp(dtype(1),'S')
        xlims = length(cc_t)./samprate - [1.5*trudata.(dtype)(1).nsamp./samprate 0];
    end
    set(xa,'fontsize',15,'xlim',xlims,'ylim',1.1*cc_max*[-1 1],'color','none')
    ylabel(xa,'\textbf{Normalised amplitude} ','fontsize',18,'interpreter','latex')
    title(xa,regexprep(dtype,'_','-'),'fontsize',22,'pos',[mean(xlims),0.85*cc_max,0])
end
fprintf('\n');


%% SW
if ~isempty(trudata.SW.phV)
axes(ax3)
plot(ax3,1./savedata.SW.periods(savedata.gdmods(idspy),:)',savedata.SW.phV(savedata.gdmods(idspy),:)','color',synthcol,'linewidth',0.2);
hp(1) = plot(ax3,1./trudata.SW.periods,trudata.SW.phV,'ok','linewidth',3,'markersize',10);
set(ax3,'fontsize',15,'xlim',[0 1./trudata.SW.periods(1)],'color','none')
xlabel(ax3,'\textbf{Frequency (Hz)}','fontsize',18,'Interpreter','latex')
ylabel(ax3,'\textbf{Phase Velocity (km/s)}','fontsize',18,'Interpreter','latex')
else
delete(ax3) 
end

%% Legend
hp(2) = plot(ax3,[-2 -1],[-2 -1],'color',synthcol,'linewidth',5);
hl = legend(ax3,hp,'Input data','trial data');
set(hl,'Fontsize',16,'Fontweight','bold')
set(hl,'pos',get(hl,'pos')+[-0.01 -0.02 +0.02 +0.04])

% 
% ps_xlims = [0 1.5*trudata.PsRF.nsamp./trudata.PsRF.samprate];
% sp_xlims = [0 1.5*trudata.SpRF.nsamp./trudata.SpRF.samprate];
% 
% 
% %% Ps
% if ~isempty(trudata.PsRF.ZRT)
%     cla(ax1)
%     for is = 1:size(savedata.PsRF.Z,2)
%         Ps_v1(:,is) = conv(trudata.PsRF(1).ZRT(:,1),savedata.PsRF.R(:,is),'full'); % Vobs*Hpre
%         Ps_v2(:,is) = conv(trudata.PsRF(1).ZRT(:,2),savedata.PsRF.Z(:,is),'full'); % Hobs*Vpre
%     end
% %     Ps_tru = conv(trudata.PsRF.ZRT(:,2),trudata.PsRF.ZRT(:,1),'full'); % Vobs*Hobs or Hobs*Vobs (identical)
%     Ps_max = max(max(abs([Ps_v1,Ps_v2])));
%     Ps_t = [0:size(Ps_v1,1)-1]'./trudata.PsRF(1).samprate;
%     plot(ax1,Ps_t,Ps_v1,'linewidth',1.5,'color',synthcol+[0 0 0.2])
%     plot(ax1,Ps_t,Ps_v2,'linewidth',1.5,'color',synthcol+[0.2 0 0]), 
% %     plot(ax1,Ps_t,Ps_tru,'linewidth',2,'color','k'), 
%     set(ax1,'fontsize',15,'xlim',ps_xlims,'ylim',1.1*Ps_max*[-1 1],'color','none')
%     ylabel(ax1,'\textbf{Normalised amplitude} ','fontsize',18,'interpreter','latex')
%     title(ax1,'Ps','fontsize',22,'pos',[mean(ps_xlims),0.80*Ps_max,0])
% end
% 
% %% Sp
% if ~isempty(trudata.SpRF.ZRT)
%     cla(ax2)
%     for is = 1:size(savedata.SpRF.Z,2)
%         Sp_v1(:,is) = conv(trudata.SpRF.ZRT(:,1),savedata.SpRF.R(:,is),'full'); % Vobs*Hpre
%         Sp_v2(:,is) = conv(trudata.SpRF.ZRT(:,2),savedata.SpRF.Z(:,is),'full'); % Hobs*Vpre
%     end
% %     Sp_tru = conv(trudata.SpRF.ZRT(:,2),-trudata.SpRF.ZRT(:,1),'full'); % Vobs*Hobs or Hobs*Vobs (identical)
%     Sp_max = max(max(abs([Sp_v1,Sp_v2])));
%     Sp_t = (size(Sp_v1,1)-[1:size(Sp_v1,1)]')./trudata.SpRF.samprate;
%     plot(ax2,Sp_t,Sp_v1,'linewidth',1.5,'color',synthcol+[0 0 0.2])
%     plot(ax2,Sp_t,Sp_v2,'linewidth',1.5,'color',synthcol+[0.2 0 0]), 
% %     plot(ax2,Sp_t,Sp_tru,'linewidth',2,'color','k'), 
%     set(ax2,'fontsize',15,'xlim',sp_xlims,'ylim',1.1*Sp_max*[-1 1],'color','none')
%     xlabel(ax2,'\textbf{Time (s)}','fontsize',18,'interpreter','latex')
%     ylabel(ax2,'\textbf{Normalised amplitude} ','fontsize',18,'interpreter','latex')
%     title(ax2,'Sp','fontsize',22,'pos',[mean(sp_xlims),0.80*Sp_max,0])
% end
% 
% 
% if iflodata
% %% Ps_lo
% if ~isempty(trudata.PsRF_lo.ZRT)
%     for is = 1:size(savedata.PsRF_lo.Z,2)
%         Ps_lo_v1(:,is) = conv(trudata.PsRF_lo.ZRT(:,1),savedata.PsRF_lo.R(:,is),'full'); % Vobs*Hpre
%         Ps_lo_v2(:,is) = conv(trudata.PsRF_lo.ZRT(:,2),savedata.PsRF_lo.Z(:,is),'full'); % Hobs*Vpre
%     end
%     Ps_lo_tru = conv(trudata.PsRF_lo.ZRT(:,2),trudata.PsRF_lo.ZRT(:,1),'full'); % Vobs*Hobs or Hobs*Vobs (identical)
%     Ps_lo_max = max(abs(Ps_lo_tru));
%     Ps_lo_t = [0:size(Ps_lo_v1,1)-1]'./trudata.PsRF_lo.samprate;
%     plot(ax4,Ps_lo_t,Ps_lo_v1,'linewidth',1.5,'color',synthcol)
%     plot(ax4,Ps_lo_t,Ps_lo_v2,'linewidth',1.5,'color',synthcol), 
%     plot(ax4,Ps_lo_t,Ps_lo_tru,'linewidth',2,'color','k'), 
%     set(ax4,'fontsize',15,'xlim',ps_xlims,'ylim',1.1*Ps_lo_max*[-1 1],'color','none')
%     title(ax4,'Ps (LOW-F)','fontsize',22,'pos',[mean(ps_xlims),0.80*Ps_lo_max,0])
% end
% 
% %% Sp_lo
% if ~isempty(trudata.SpRF_lo.ZRT)
%     for is = 1:size(savedata.SpRF_lo.Z,2)
%         Sp_lo_v1(:,is) = conv(trudata.SpRF_lo.ZRT(:,1),savedata.SpRF_lo.R(:,is),'full'); % Vobs*Hpre
%         Sp_lo_v2(:,is) = conv(trudata.SpRF_lo.ZRT(:,2),savedata.SpRF_lo.Z(:,is),'full'); % Hobs*Vpre
%     end
%     Sp_lo_tru = conv(trudata.SpRF_lo.ZRT(:,2),trudata.SpRF_lo.ZRT(:,1),'full'); % Vobs*Hobs or Hobs*Vobs (identical)
%     Sp_lo_max = max(abs(Sp_lo_tru));
%     Sp_lo_t = (length(Sp_lo_v1)-[1:length(Sp_lo_v1)]')./trudata.SpRF_lo.samprate;
%     plot(ax5,Sp_lo_t,Sp_lo_v1,'linewidth',1.5,'color',synthcol)
%     plot(ax5,Sp_lo_t,Sp_lo_v2,'linewidth',1.5,'color',synthcol), 
%     plot(ax5,Sp_lo_t,Sp_lo_tru,'linewidth',2,'color','k'), 
%     set(ax5,'fontsize',15,'xlim',sp_xlims,'ylim',1.1*Sp_lo_max*[-1 1],'color','none')
%     xlabel(ax5,'\textbf{Time (s)}','fontsize',18,'interpreter','latex')
%     title(ax5,'Sp (LOW-F)','fontsize',22,'pos',[mean(sp_xlims),0.80*Sp_lo_max,0])
% end
% 
% end


pause(0.001)


if ifsave
    fprintf('Saving...')
    save2pdf(61,ofile,'/');
    fprintf('\n')
end

end


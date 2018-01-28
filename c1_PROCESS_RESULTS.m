function [misfits,allmodels,goodchains] = c1_PROCESS_RESULTS( misfits,allmodels,par,ifsave,ofile )
% [misfits,allmodels,goodchains] = c1_PROCESS_RESULTS( misfits,allmodels,par,ifsave,ofile )
% 
% Script to process the results and make some plots of the misfit + the
% likelihood changing with iteration

if nargin < 4 || isempty(ifsave)
    ifsave=false;
end
if nargin < 5 || isempty(ofile)
    ofile = 'figs/misfits_vs_iter';
end


figure(88);clf, set(gcf,'pos',[850 198 900 888])
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2); hold(ax2,'on');

pdm_min = inf;
pdm_max = -inf;
ch2_min = inf;
ch2_max = -inf;

downsampfac = 2;



%% loop through each chain
nchains = length(misfits);
goodchains = true(nchains,1);

for iii = 1:nchains
    
    if nchains>1, mf = misfits{iii};   else mf = misfits; end
    if nchains>1, am = allmodels{iii}; else am = allmodels; end
    if isempty(am), continue; end
    basecol = colour_get(iii,nchains+1,0,parula); basecol = basecol(:)';
    
    %% TRIM MISFITS AND ALLMODELS STRUCTURES

    N = am(1).Nstored;

    % kill allmodels beyond N
    am(N+1:end) = [];
    

    % kill misfits beyond N
    fns = fieldnames(mf);
    for ifn = 1:length(fns)
        if length(mf.(fns{ifn})) > N
            mf.(fns{ifn})(N+1:end) = [];
        end
    end
    
    %% move on if not enough good results
    if length(am) < par.inv.bestNmod2keep && ~isinf(par.inv.bestNmod2keep)
        mf.bestmods = false(mf.Nstored,1);
        am = dealto(am,'bestmods',false);
        goodchains(iii) = false;
        if nchains>1, misfits{iii} = mf;   else misfits = mf; end
        if nchains>1, allmodels{iii} = am; else allmodels = am; end
        continue
    end
    
    %% identify "best" models 
%     from those in the top bestNmod2keep based on summed ranking of their fits to each datatype

    bestind = [1:length(mf.iter)]';
    bestind(mf.iter <= par.inv.burnin) = []; % kill all before burnin

    if ~isinf(par.inv.bestNmod2keep) % subset if not inf to keep. Else keep all

    if par.inv.bestNmod2keep>0 % if specifying how many to keep based on low error

        score_overall = zeros(length(bestind),length(par.inv.datatypes));
        for id = 1:length(par.inv.datatypes);
            dtype = par.inv.datatypes{id};
            pdtyp = parse_dtype(dtype);
    %         if strcmp(pdtype{1},'BW') && (~strcmp(pdtype{3},'def') || ~strcmp(pdtype{4},'def')), continue; end
            chi2 = [mf.chi2.(dtype)]';
            [~,irank_mf] = sort(sum(chi2(bestind,:),2));
            [~,score_overall(:,id)] = sort(irank_mf);
        end
        score_overall = sum(score_overall,2);

    %     if isfield(mf.chi2,'BW_Ps') && isfield(mf,'BW_Ps') && isfield(mf,'SW')
    %         [~,irank_mf_ps] = sort(mf.chi2(bestind)); [~,rank_mf_ps] = sort(irank_mf_ps);
    %         [~,irank_mf_sp] = sort(mf.chi2_sp(bestind)); [~,rank_mf_sp] = sort(irank_mf_sp);
    %         [~,irank_mf_sw] = sort(mf.chi2_SW(bestind)); [~,rank_mf_sw] = sort(irank_mf_sw);
    %         score_overall = rank_mf_ps + rank_mf_sp + rank_mf_sw;
    %     else
    %             [~,irank_mf] = sort(mf.chi2sum(bestind)); [~,score_overall] = sort(irank_mf);
    %     end

        sort_score_overall = sort(score_overall);
        min_score_overall = sort_score_overall(min([par.inv.bestNmod2keep,length(bestind)]));
        bestind(score_overall>min_score_overall) = [];
    elseif par.inv.bestNmod2keep<0
        bestind = bestind(randperm(length(bestind),min([-par.inv.bestNmod2keep,length(bestind)])));
    end
    end
    
    mf.bestmods = false(mf.Nstored,1);
    mf.bestmods(bestind) = true;
    
%     subplot(211),semilogy(mf.iter,mf.chi2,'bo',mf.iter(mf.bestmods),mf.chi2(mf.bestmods),'ro')
%     bestind = [1:mf.Nstored]';
%     bestind(misfits.iter <= par.inv.burnin) = [];
%     [~,rank_mf] = sort(mf.chi2(bestind));
%     max_rank_ind = min([par.inv.bestNmod2keep,length(bestind)]);
%     bestind = bestind(rank_mf(1:max_rank_ind));
%     mf.bestmods = false(mf.Nstored,1);
%     mf.bestmods(bestind) = true;
%     subplot(212),semilogy(mf.iter,mf.chi2,'bo',mf.iter(mf.bestmods),mf.chi2(mf.bestmods),'ro')

    am = dealto(am,'bestmods',mf.bestmods);

    %% plug back into structures
    if nchains>1, misfits{iii} = mf;   else misfits = mf; end
    if nchains>1, allmodels{iii} = am; else allmodels = am; end

    
    %% PLOT IMPROVEMENT IN MODEL FIT
    downsamp = (1:downsampfac:length(am));

    [ax1,h1,h2] = plotyy(ax1,mf.iter(downsamp),mf.chi2sum(downsamp),mf.iter(downsamp),mf.logLike(downsamp),'semilogy');
    
    hold(ax1(1),'on'),hold(ax1(2),'on')
    pdm_min = min([pdm_min,min(mf.logLike)-1]);
    pdm_max = max([pdm_max,max(mf.logLike)+10]);
    ch2_min = min([ch2_min,min(mf.chi2sum)/2]);
    ch2_max = max([ch2_max,max(mf.chi2sum)*2]);
%     set(ax1(1),'Yscale','log','ylim',[min(mf.chi2)/2 max(mf.chi2)*2],'ytick',round_level(linspace(min(mf.chi2)/2,max(mf.chi2)*2,5),5));
%     ytk = unique(round_level(linspace(min(mf.logLike),max(mf.logLike),5),5));
%     set(ax1(2),'Yscale','log','ylim',[min(mf.logLike)-1 max(mf.logLike)+10],...
%         'ytick',ytk,'yticklabel',ytk);
    set(ax1(1),'Yscale','log','ylim',[ch2_min ch2_max],'ytick',round_level(linspace(ch2_min,ch2_max,5),5));
    ytk = unique(round_level(linspace(pdm_min+1,pdm_max-10,5),5));
    set(ax1(2),'Yscale','log','ylim',[pdm_min pdm_max],...
        'ytick',ytk,'yticklabel',ytk);
    
    set(h1,'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',[0 0.447 0.741]);
    set(h2,'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',[0.85 0.325 0.098]);
    set(get(ax1(1),'ylabel'),'String','$\chi^2$ misfit','Fontsize',20,'interpreter','latex')
    set(get(ax1(2),'ylabel'),'String','$\log_{10}{\,p(m|d)}$','Fontsize',20,'interpreter','latex')
    set(ax1,'Fontsize',16)
    
    
    cols = [[0 0.447 0.741];[0.85 0.325 0.098];[0.929 0.694 0.125];[0.494 0.184 0.556];[0.466 0.674 0.188];[0.1 0.76 0.288]];
    
    for idt = 1:length(par.inv.datatypes) 
        dtype = par.inv.datatypes{idt};
        rms = [mf.rms.(dtype)]';
        hrms(idt)=plot(ax2,mf.iter(downsamp),sum(rms(downsamp,:),2),'o'); hold on
        set(hrms(idt),'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',cols(idt,:));
    end
    hrms(idt+1) = hrms(1); set(hrms(idt+1),'markerfacecolor','none');
    set(hrms(2:idt+1),'linewidth',2)
    set(ax2,'yscale','log','Fontsize',16 )
    ylabel(ax2,'RMS misfit','Fontsize',20,'interpreter','latex')
    xlabel(ax2,'Iteration','Fontsize',20)
    hl = legend(hrms(1:idt),strrep(par.inv.datatypes,'_','-'));

    
    % rate of acceptance past burnin

end

%% title
htit = title_custom([par.sta,' ',par.nwk],0.5,'fontweight','bold','fontsize',25);


%% determine "good" models

if ifsave
    fprintf('saving, may take a while\n')
    save2png(88,ofile,'/');
end
pause(0.05)

%% PARSE GOOD AND BAD CHAINS
rms_alldata = nan(nchains,length(par.inv.datatypes));

% gather average rms errors of each chain
for iii=1:nchains
    if goodchains(iii)==false, rms_alldata(iii,:) = nan; continue; end % already know it's bad
    for id = 1:length(par.inv.datatypes)
        % assign structures
        if nchains>1, mf = misfits{iii}; else mf = misfits; end
        if isempty(mf), goodchains(iii)=false; continue; end
        ind = mf.iter > par.inv.burnin;
        dtype = par.inv.datatypes{id};
        % work out rms for this dtype (average across all data streams for this dtype)
        rms = [mf.rms.(dtype)]';
        rms_alldata(iii,id) = mean(sum(rms(ind,:),2)); 
        % work out if the chain got stuck - if there is no change to the
        % data over many iterations - must be stuck for 500 iterations to
        % signify
        if any(any(diff(rms(ind,:),ceil(1000./par.inv.saveperN),1)==0))
            rms_alldata(iii,id) = nan; continue;
        end
    end
end

% goodchains = true(nchains,1);
for id = 1:length(par.inv.datatypes)
    goodchains = goodchains & (rms_alldata(:,id) < 1.3*nanmean(rms_alldata(:,id)));
end
goodchains=find(goodchains);

%% GET TARGET MODEL for comparison
global TRUEmodel
if ~isempty(TRUEmodel)
Z = TRUEmodel.Z;
vs = TRUEmodel.vs;
vp = TRUEmodel.vp;
rho = TRUEmodel.rho;
fprintf('TRUE sed thickness = %.1f km\n',Z(find(vs>=3.2,1,'first')))
fprintf('TRUE moho depth = %.1f km\n',Z(find(vs>4.0,1,'first')))
fprintf('TRUE Vs seds top = %.1f km/s\n',vs(1))
fprintf('TRUE Vs seds bot = %.1f km/s\n',vs(find(vs<3.45,1,'last')))
fprintf('TRUE Vs crust top = %.1f km/s\n',vs(find(vs>3.45,1,'first')))
fprintf('TRUE Vs crust bot = %.1f km/s\n',vs(find(rho<3.2,1,'last')))
fprintf('TRUE fractional dVs sed/crust = %.1f %% \n',-100*(vs(find(vs<3.5,1,'last')) - vs(find(vs>3.5,1,'first')))/vs(find(vs<3.5,1,'last')))
fprintf('TRUE fractional dVs crust/mantle = %.1f %% \n',-100*(vs(find(rho<3.2,1,'last')) - vs(find(rho>3.2,1,'first')))/vs(find(rho<3.2,1,'last')))
for ii = linspace(par.mod.sed.hmax+par.mod.crust.hmax,par.mod.maxz,6)
    try fprintf('TRUE Vs at %.0f km = %.2f km/s\n',ii,linterp(Z,vs,ii));end
end
end

% 
% %% PLOT PROGRESS OF MODEL
% figure(99); clf; set(gcf,'pos',[120 151 920 947])
% Nacc = allmodels.Nstored;
% for ii = 1:3:Nacc 
%     if allmodels(ii).iter <= par.inv.burnin, continue, end
%     
%     if mod(ii,30)==1, 
%     figure(99);
%     subplot(131), hold on; 
%     plot(vs,Z,'b','Linewidth',2)
%     subplot(132), hold on;
%     plot(vp,Z,'b','Linewidth',2)
%     subplot(133), hold on;
%     plot(rho,Z,'b','Linewidth',2)
%     end
%     
%     subplot(131), hold on;
%     hp1 = plot(allmodels(ii).VS,allmodels(ii).z,'-o',...
%         'color',colour_get(misfits.chi2(ii),50,1,autumn));
%     set(gca,'ydir','reverse');
%     
%     subplot(132), hold on;
%     hp2 = plot(allmodels(ii).VP,allmodels(ii).z,'-o',...
%         'color',colour_get(misfits.chi2(ii),50,1,autumn));
%     set(gca,'ydir','reverse')
% 
%     subplot(133), hold on;
%     hp3 = plot(allmodels(ii).rho,allmodels(ii).z,'-o',...
%         'color',colour_get(misfits.chi2(ii),50,0,autumn));
%     set(gca,'ydir','reverse')
%     
%     
%     pause(0.01)
%     figure(99);
%     set(hp1,'color',[0.5 0.5 0.5]),set(hp2,'color',[0.5 0.5 0.5]),set(hp3,'color',[0.5 0.5 0.5])
% end


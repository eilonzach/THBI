function [ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels,ofile )
% [ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels,ofile )

if nargin<2
    ofile=[];
end

downsampfac = 15;

figure(56), clf, set(gcf,'pos',[65 682 814 423]), hold on

% dtypcols = [[0 0.447 0.741];[0.85 0.325 0.098];[0.929 0.694 0.125];[0.494 0.184 0.556];[0.466 0.674 0.188];[0.1 0.76 0.288]];
    
    strmmark = {'o','v','^','s','p'};

%% loop through each chain
if iscell(allmodels),nchains = length(allmodels); else nchains=1; end
for iii = 1:nchains
    basecol = colour_get(iii,nchains+1,0,parula); basecol = basecol(:)';
    
    if nchains>1
        am = allmodels{iii};
        if isempty(am); continue; end
    else 
        am = allmodels; 
        if iscell(am), am = am{1}; end
    end

    dp = [am.datahparm]; 
    
    downsamp = (1:downsampfac:length(am));

    dtypes = fieldnames(dp);
    Nd = length(dtypes);
    dtypcols = colour_get([1:Nd],Nd,1,[parula;flipud(spring)]);

    clear('hypparm_trends');
    for id = 1:length(dtypes)
        Nt = length(dp(1).(dtypes{id}));    
        hpt = reshape([dp.(dtypes{id})],Nt,length(dp))';
        for it = 1:Nt
            hypparm_trends(:,id,it) = hpt(:,it);

            hp(id) = plot([am(downsamp).iter],hpt(downsamp,it));

            set(hp(id),'marker',strmmark{it},'linestyle','none','markersize',5,...
                'markerfacecolor',basecol,'markeredgecolor',dtypcols(id,:));
        end
    end
    
    set(gca,'yscale','log','fontsize',16)
    xlabel('Iteration','fontsize',20,'interpreter','latex')
    ylabel('$\log(\sigma_i)$','fontsize',20,'interpreter','latex')
    title('Hyperparameter variation','fontsize',20)
    dlabs = regexprep(regexprep(dtypes,'_','-'),'sigma','');
    hl = legend(hp,dlabs,'location','northeast','fontsize',17);
    
end 

if ~isempty(ofile)
    save2pdf(56,ofile,'/')
end


end
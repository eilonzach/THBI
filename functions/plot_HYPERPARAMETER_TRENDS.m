function [ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels )
% [ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels )

figure(56), clf, set(gcf,'pos',[65 682 814 423]), hold on

%% loop through each chain
if iscell(allmodels),nchains = length(allmodels); else nchains=1; end
for iii = 1:nchains
    basecol = colour_get(iii,nchains+1,0,parula); basecol = basecol(:)';
    
    if nchains>1, 
        am = allmodels{iii};
    else 
        am = allmodels;
    end

    dp = [am.datahparm]; 

    dtypes = fieldnames(dp);
    clear('hypparm_trends');
    for id = 1:length(dtypes)
        hypparm_trends(:,id) = [dp.(dtypes{id})]';
    end

  
    h1 = plot([am.iter],hypparm_trends(:,1),'o');
    h2 = plot([am.iter],hypparm_trends(:,2),'o');
    h3 = plot([am.iter],hypparm_trends(:,3),'o');
    h4 = plot([am.iter],hypparm_trends(:,4),'o');
    h5 = plot([am.iter],hypparm_trends(:,5),'o');

    set(h1,'marker','o','linestyle','none','markersize',5,...
       'markerfacecolor',basecol,'markeredgecolor',[0 0.447 0.741]);
    set(h2,'marker','o','linestyle','none','markersize',5,...
       'markerfacecolor',basecol,'markeredgecolor',[0.85 0.325 0.098]);
    set(h3,'marker','o','linestyle','none','markersize',5,...
       'markerfacecolor',basecol,'markeredgecolor',[0.929 0.694 0.125]);
    set(h4,'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',[0.494 0.184 0.556]);
    set(h5,'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',[0.466 0.674 0.188]);
    set(gca,'yscale','log','fontsize',16)
    xlabel('Iteration','fontsize',20,'interpreter','latex')
    ylabel('$\log(\sigma_i)$','fontsize',20,'interpreter','latex')
    title('Hyperparameter variation','fontsize',20)
    dlabs = regexprep(regexprep(dtypes,'_','-'),'sigma','');
    hl = legend([h1,h2,h3,h4,h5],dlabs,'location','northeast','fontsize',17);
    
end 

end
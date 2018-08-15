function plot_quickmodel(par,model,model1)
% plot_quickmodel(par,model,model1)

global TRUEmodel

if isempty(gcp('nocreate')) && par.inv.verbose 
        
    figure(85);clf,set(gcf,'pos',[1194 4 622 529])
    ax1 = subplot(1,2,1); hold on, 
    ax2 = subplot(1,2,2); hold on, 

    if ~isempty(TRUEmodel)
    plot(ax1,TRUEmodel.VS,TRUEmodel.z,'k','linewidth',1); 
    plot(ax2,TRUEmodel.VP,TRUEmodel.z,'k','linewidth',1);
    end
    
    plot(ax1,model.VS,model.z,'r','linewidth',1.5);
    plot(ax2,model.VP,model.z,'r','linewidth',1.5);

    plot(ax1,model1.VS,model1.z,'b--','linewidth',1.5);    
    plot(ax2,model1.VP,model1.z,'b--','linewidth',1.5);
    
    set(ax1,'ydir','reverse');
    set(ax2,'ydir','reverse');
    pause(0.01)
        
end

end





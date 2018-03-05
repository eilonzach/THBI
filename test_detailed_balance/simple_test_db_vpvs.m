Niter=5e6; 
a = nan(Niter,1); b = nan(Niter,1); 
a(1) = random('unif',1.6,1.9,1); 
for ii=2:Niter
    da = random('norm',0,0.01); 
    if (a(ii-1)+da)>=1.6 && (a(ii-1)+da)<=1.9
        a(ii) = a(ii-1)+da; b(ii) = 1;
    else
        a(ii) = a(ii-1); b(ii) = 0;
    end
end

figure(44);
X = midpts(linspace(1.6,1.9,40));
N1 = hist(a,X);
N2 = hist(a(b==1),X);
subplot(211), hold on, 
bar(X,N1./sum(N1),'facecolor',[0.9 0.1 0.1],'edgecolor','k','BarWidth',1),xlim([1.6,1.9]) 
plot([1.6,1.9],[1,1]./length(X),'--b','linewidth',2)
subplot(212), hold on
bar(X,N2./sum(N2),'facecolor',[0.9 0.1 0.1],'edgecolor','k','BarWidth',1), xlim([1.6,1.9]) 
plot([1.6,1.9],[1,1]./length(X),'--b','linewidth',2)


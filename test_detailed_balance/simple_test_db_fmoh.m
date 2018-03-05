clear
Niter=1e6; 
c = nan(Niter,1); m = nan(Niter,1); 
dv = 50;

while dv>30 || dv <= 0
    c(1) = random('unif',2.5,4.3,1); 
    m(1) = random('unif',3.7,5.1,1); 
    dv = 200*(m(1)-c(1))./(m(1)+c(1));
end
for ii=2:Niter
    dd = random('norm',0,0.08); 
    id = randi(2);
    
    if     id==1, c_ = c(ii-1)+dd; m_ = m(ii-1);
    elseif id==2, c_ = c(ii-1); m_ = m(ii-1)+dd;
    end
    
    dv = 200*(m_-c_)./(m_+c_);
    
    if c_>=2.5 && c_<=4.3 && m_>=3.7 && m_<=5.1 && dv<=30 && dv>0
        c(ii) = c_; 
        m(ii) = m_;
    else
        c(ii) = c(ii-1); 
        m(ii) = m(ii-1);
    end
end

figure(44);clf
Xc = midpts(linspace(2.5,4.3,20));
Xm = midpts(linspace(3.7,5.1,20));
Xf = midpts(linspace(0,30,20));
N1 = hist(c,Xc);
N2 = hist(m,Xm);
N3 = hist(200*(m-c)./(m+c),Xf);
subplot(311), hold on, 
bar(Xc,N1./sum(N1),'facecolor',[0.9 0.1 0.1],'edgecolor','k','BarWidth',1),xlim([2.5,4.3]) 
subplot(312), hold on
bar(Xm,N2./sum(N2),'facecolor',[0.9 0.1 0.1],'edgecolor','k','BarWidth',1), xlim([3.7,5.1]) 
subplot(313), hold on
bar(Xf,N3./sum(N3),'facecolor',[0.9 0.1 0.1],'edgecolor','k','BarWidth',1), xlim([0,30]) 


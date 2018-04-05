Niter=2e5; 
a = nan(Niter,1); 
Nmax = 15;
Nmin = 3;
a(1) = random('unid',Nmax-Nmin+1,1) + Nmin - 1 ; 

kk=1;
for ii=2:Niter, bd=true; 
    Nnew = a(kk) + (-1)^random('unid',2,1); % add or subtract
    p_bd = a(kk)/Nnew;
    if Nnew<Nmin || Nnew>Nmax
        kk = kk+1;
        a(kk) = a(kk-1);
        continue
    end
    if p_bd < rand
        kk = kk+1;
        a(kk) = a(kk-1);
        continue
    end
    % if at this point, within bounds and accepted
    kk = kk+1;
    a(kk) = Nnew;
end

X = [Nmin:1:Nmax];
N = hist(a,X)/length(a);
predN = (1./X)./sum(1./X);

figure(44); clf; hold on

bar(X,N),xlim([Nmin-0.5 Nmax+0.5]);
plot(X,predN,'-r')
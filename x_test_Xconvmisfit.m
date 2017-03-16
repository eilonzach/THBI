% PS MISFIT
Vobs = trudata.PsRF.ZRT(:,1);
Hobs = trudata.PsRF.ZRT(:,2);
Vpre = predata.PsRF.ZRT(:,1);
Hpre = predata.PsRF.ZRT(:,2);

% Vpre = trudata.PsRF.ZRT(:,1);
% Hpre = trudata.PsRF.ZRT(:,2);

sumEobs = Vobs'*Vobs + Hobs'*Hobs;
sumEpre = Vpre'*Vpre + Hpre'*Hpre;

%% Standard xconv method

e1 = conv(Vobs,Hpre,'full') - conv(Hobs,Vpre,'full');
% e = conv(Vobs,Hpre,'same') - conv(Hobs,Vpre,'same');

misfit1 = e1'*e1

%% Wraparound + xconv 
nh = length(Hpre);
cv1 = conv([Vobs;Vobs],Hpre,'full');
cv1 = cv1(nh:2*nh);
cv2 = conv([Hobs;Hobs],Vpre,'full');
cv2 = cv2(nh:2*nh);

e2 = cv1 - cv2;
misfit2 = e2'*e2

%% fourier domain
[ fHpre,ffhp ] = fft_ze( Hpre );
[ fVpre,ffvp ] = fft_ze( Vpre );
[ fHobs,ffho ] = fft_ze( Hobs );
[ fVobs,ffvo ] = fft_ze( Vobs );

fcv3 = fVobs.*fHpre;
cv3 = ifft(fcv3);
fcv4 = fHobs.*fVpre;
cv4 = ifft(fcv4);

e3 = cv3 - cv4;
misfit3 = e3'*e3

% for k = 1:length(Hpre)
%     dhp(2*k-1) = Hpre(k);
%     dhp(2*k) = 0;
%     dvp(2*k-1) = Vpre(k);
%     dvp(2*k) = 0;    
% end
% for k = 1:length(Hobs)
%     dho(2*k-1) = Hobs(k);
%     dho(2*k) = 0;
%     dvo(2*k-1) = Vobs(k);
%     dvo(2*k) = 0;    
% end
% [ fdhp,ffdhp ] = fft_ze( dhp );
% 
% for k = 1:floor(length(Hpre)/2) + 1
%     fdhp2(k) = fdhp(2*k-1);
%     ffdhp2(k) = 2*ffdhp(2*k-1);
% end
% 
% 
% 
% 
% figure(4), clf, hold on
% Nf = floor(length(Hpre)/2) + 1;
% plot(ff(1:Nf),abs(fHpre(1:Nf)))
% % plot(ffdhp(1:Nf),abs(fdhp(1:Nf)))
% plot(ffdhp2(1:Nf),abs(fdhp2(1:Nf)))
% 

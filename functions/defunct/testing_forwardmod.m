% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with two traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
% clear all
close all
addpath('matguts')
% parms
samprate = 100;
%% trace 1 Q & V - faster c0 and less attenuated
Q0_1 = 150;
c0_1 = 4.00e3; % reference velocity in km/s
%% trace 2 Q & V - slower c0 and more attenuated
Q0_2 = 30;
c0_2 = 3.95e3; % reference velocity in km/ss

L = 200e3;

alpha = 0;

% calc. true values
dT = L*(1./c0_2 - 1./c0_1);
dtstar = L*(1./c0_2./Q0_2 - 1./c0_1./Q0_1);
fprintf('Theoretical dT=%.2f and dtstar=%.2f\n',dT,dtstar)

%% ========== Make the two traces ==========

%% Make simple pulse
dt = 1/samprate;
tt = [-50:dt:50]';
T = tt(end)-tt(1);
N = T/dt; tt = tt(1:N);
fnq = 0.5/dt;
dat = synthtrace(T,2,1,dt,'gauss');

%% Delay pulses
dT1 = L*(1./c0_1 - 1./mean([c0_1,c0_2]));
dat1 = interp1(tt+dT1,dat,tt,'linear',0)';
dT2 = L*(1./c0_2 - 1./mean([c0_1,c0_2]));
dat2 = interp1(tt+dT2,dat,tt,'linear',0)';

%% Attenuate pulses
% take fft of pulse
[DAT_1,~ ] = fft_ze(dat1,dt);
[DAT_2,ff] = fft_ze(dat2,dt);
w = 2*pi*ff;
% Make attenuation operators
[ Dwt_1 ] = attenuation_operator( Q0_1,c0_1,L,w,alpha);
[ Dwt_2 ] = attenuation_operator( Q0_2,c0_2,L,w,alpha);
% apply attenuation operator and ifft
qdat1 = abs(ifft(DAT_1.*Dwt_1));
qdat2 = abs(ifft(DAT_2.*Dwt_2));

%% Plot original and attenuated pulses
figure(1), clf, hold on
plot(tt,dat,'k','LineWidth',1.5)
plot(tt,qdat1,'b','LineWidth',1.5)
plot(tt,qdat2,'r','LineWidth',1.5)
title('Original (black), trace1 (blue) and trace2 (red)')

%% Attenuate 1 to try to get 2

odat1 = attenuate( qdat1,samprate,dtstar,dT,1.1,0 );

% plot
figure(1), hold on
plot(tt,odat1,'--g','LineWidth',1.5)

% calc misfit
misfit = sum((odat1-qdat2).^2);

%% Make misfit map
dtstar_test = [0:0.02:2]';
dT_test = [0:0.02:1]';

M1 = length(dtstar_test);
M2 = length(dT_test);


E = zeros(M1,M2);
for ii = 1:M1
for jj = 1:M2
    odat1 = attenuate( qdat1,samprate,dtstar_test(ii),dT_test(jj),1,0 );
    E(ii,jj) = sum((odat1-qdat2).^2);
end
end

[Emin,xi,yi] = mingrid(E);
dT_est = dT_test(xi);
dtstar_est = dtstar_test(yi);

figure(2), clf, hold on
pcolor(dT_test,dtstar_test,E)
plot(dT,dtstar,'ow','Linewidth',2)
plot(dT_est,dtstar_est,'or','Linewidth',2)





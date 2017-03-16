function [data,tt,samprate,dT_tru,dtstar_tru] = testdata1(ifplot)
% [data,samprate,dT_tru,dtstar_tru] = testdata1(ifplot)
% 
% test dataset 1

if nargin==0
    ifplot=1;
end



% parms
samprate = 40;
T = 1e3;

%% trace 1 Q & V - faster c0 and less attenuated
Q0_1 = 1000;
c0_1 = 4.0e3; % reference velocity in km/s
%% trace 2 Q & V - slower c0 and more attenuated
Q0_2 = 150;
c0_2 = 3.95e3; % reference velocity in km/s
%% trace 3 Q & V - faster c0 and less attenuated
Q0_3 = 50;
c0_3 = 3.92e3; % reference velocity in km/s

L = 100e3;

alpha = 0.09;


% calc. true values
TT = L*[c0_1,c0_2,c0_3].^-1;
dT_tru = demean(TT)';

tstar = L*[c0_1*Q0_1,c0_2*Q0_2,c0_3*Q0_3].^-1;
dtstar_tru = demean(tstar)';


%% ========== Make the two traces ==========

%% Make simple pulse
dt = 1/samprate;
tt = [-T/2:dt:T/2]'; % well padded
T = tt(end)-tt(1);
N = T/dt; tt = tt(1:N);
fnq = 0.5/dt;
dat = synthtrace(T,4,1,dt,'gauss');


%% Attenuate pulses
% take fft of pulse
[DAT,ff] = fft_ze(dat',dt);

w = 2*pi*ff;

% Make attenuation operators
[ Dwt_1 ] = attenuation_operator( Q0_1,c0_1,L,w,alpha,'zph');
[ Dwt_2 ] = attenuation_operator( Q0_2,c0_2,L,w,alpha,'zph');
[ Dwt_3 ] = attenuation_operator( Q0_3,c0_3,L,w,alpha,'zph');

% apply attenuation operator and ifft
qdat1 = real(ifft(DAT.*Dwt_1));
qdat2 = real(ifft(DAT.*Dwt_2));
qdat3 = real(ifft(DAT.*Dwt_3));

%% Delay pulses
qdat1 = interp1(tt+dT_tru(1),qdat1,tt,'linear',0); qdat1 = qdat1(:);
qdat2 = interp1(tt+dT_tru(2),qdat2,tt,'linear',0); qdat2 = qdat2(:);
qdat3 = interp1(tt+dT_tru(3),qdat3,tt,'linear',0); qdat3 = qdat3(:);

qdat = [qdat1,qdat2,qdat3];

if ifplot
%% Plot original and attenuated pulses
figure(1), clf, hold on
plot(tt,dat,'k','LineWidth',1.5)
plot(tt,qdat,'LineWidth',1.5)
title('Original (black), and five delayed traces')
xlim([-10 10])
end

%% output
data=qdat;
nstas = size(data,2);


end

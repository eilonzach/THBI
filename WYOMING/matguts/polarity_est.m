function [ polests,polests_cc,polests_maxab,polests_mami ] = polarity_est( dat_template,dat,maxlag )
%[ polests ] = polarity_est( dat_template,dat )
%   Function to extimate the polarity of a bunch of arrivals in time series
%   that comprise columns of the NptxNtr matrix "dat". Each of these time
%   series is, in turn, compared to "dat_template" using cross-correlation.
%   "dat_template" is some idealised time series - either use a simple
%   gaussian pulse, or use the highest SNR trace in a record section. 
%   
%   The function then cross correlates dat_template with each column of dat
%   and then with that same column but with flipped polarity. The polarity
%   that gives the highest max cross correlation is then output. "maxlag"
%   (in units of # samples, lag_seconds*samplerate) defines the limits of
%   the cross correlation lags. This gives us one estimate of the polarity:
%   polests_cc. However, sometimes (because of sinusoidal waveforms) the
%   opposite polarity is nearly as high a cc value, so only use this
%   estimate if it is well above the rms cc. Otherwise prefer using the
%   maximum absolute amplitude.
%   If the two estimates disagree, return "0" - don't use this one. 


Ntr = size(dat,2);

polests = zeros(Ntr,1);
polests_cc = zeros(Ntr,1);
polests_maxab = zeros(Ntr,1);
polests_mami = zeros(Ntr,1);

for ii = 1:Ntr
    %% measure
    cc = xcorr(dat_template,dat(:,ii),maxlag);
    polests_cc(ii) = unique(sign(maxab(cc)));
    polests_maxab(ii) = sign(maxab(dat(:,ii)))./sign(maxab(dat_template));
    [mmo,sig_mami] = min_max_ord(dat(:,ii));
    polests_mami(ii) = mmo./min_max_ord(dat_template);
    
    %% significances:
    % only use the cc estimate if the cc peak is 'significant'
%     cc_snr = abs(maxab(cc))/rms(cc);
    sig_cc = (abs(max(cc))-abs(min(cc)))./rms(cc);
    if abs(sig_cc) < 0.2, polests_cc(ii) = 0; end
    if any(sig_mami<4), polests_mami(ii)=0; end
    
    %% plotting
%     figure(2), clf
%     subplot(211), plot([dat_template/max(abs(dat_template)),...
%                         dat(:,ii)/max(abs(dat(:,ii)))])
%     subplot(212), plot(cc)
%     xlabel(['cc: ',num2str(polests_cc(ii)),'  maxab: ',num2str(polests_maxab(ii)),'  mamio: ',num2str(polests_mami(ii))],'fontsize',22)

end

polests = polests_cc + polests_maxab +polests_mami;
polests(abs(polests)<2) = 0; % -1/1 if |score| is ?2, 0 otherwise
polests = sign(polests);

end

function [mmo,sig] = min_max_ord(xx)
    [~,imi] = min(xx);
    [~,ima] = max(xx);
    mmo = sign(imi-ima); % negative if minimum comes first, positive if max comes first)
%     sig = (max(xx)-min(xx))/rms(xx);
    sig = [max(xx),-min(xx)]/rms(xx);
end
    


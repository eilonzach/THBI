function [ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par, noise_sta_deets )
%[ noisy_trudata,par ] = z2_NOISIFY_SYNTH( trudata, par, noise_sta_deets )
% 
% Function to add realistic noise to the synthetic data by grabbing an
% excerpt of noisy data from before the body wave arrivals (for Ps, Sp)
% from a real station, computing the power spectra, and adding noise to the
% synthetic data with this power spectrum (but randomised phase). 
% 
% This should be done before any of the data cleaning steps, as would be
% done for the real data. 

nwin = [-80 -30];

% correct station name if needed
sy = regexp(noise_sta_deets.sta,'SYNTH_');
if ~isempty(sy)
    noise_sta_deets.sta = noise_sta_deets.sta(sy+6:end);
end

% %% Get real data
% wd = pwd;
% cd(noise_sta_deets.datadir);
% datfiles = {};
% for ig = 1:length(noise_sta_deets.gc)
% for gcwind = -4:4;
%     df = dir(sprintf('avar_dat_%s_%s_%.0f_*.mat',noise_sta_deets.sta,noise_sta_deets.nwk,noise_sta_deets.gc(ig)+gcwind));
%     datfiles = [datfiles;{df.name}']; %#ok<AGROW>
% end
% end
% for id = 1:length(datfiles)
% avar = load(datfiles{id});
% avars(id) = avar.avar;
% end
% cd(wd);

%% assemble real data into a 


dtypes = fieldnames(trudata);
for id = 1:length(dtypes)
    dtype = dtypes{id};
    pdtyp = parse_dtype(dtype);
    if strcmp(pdtyp{1},'SW'), continue; end
    
    for ia = 1:length(trudata.(dtype))

        % grab real data
        wd = pwd;
        cd(noise_sta_deets.datadir);
        datfiles = {};
        for gcwind = -4:4
            df = dir(sprintf('avar_dat_%s_%s_%.0f_*.mat',noise_sta_deets.sta,noise_sta_deets.nwk,trudata.(dtype)(ia).gcarc+gcwind));
            datfiles = [datfiles;{df.name}']; %#ok<AGROW>
        end
        load(datfiles{1});
        cd(wd);

        ip = find(strcmp(pdtyp{2},avar.phases));

        nind = avar.tt(:,ip)>=nwin(1) & avar.tt(:,ip)<nwin(2);
        ndt = diff(avar.tt(1:2));
        idt = 1./trudata.(dtype)(1).samprate;
    
        
        odat = zeros(size(trudata.(dtype)(ia).PSV));
        for ic = 1:2
        idat = trudata.(dtype)(ia).PSV(:,ic);
        nsamps = length(idat);
        %% process noise series
        ndat = avar.dataPSVSH(nind,ic,ip);
        ndat = ndat(~isnan(ndat));
        ndat=detrend(ndat,'constant'); %detrend
%         ndat=ndat./sqrt(mean(ndat.^2)); %denorm
        %% spectrum
        nw=4;
        freq=[1:floor(nsamps/2)+1]'./(nsamps*idt);%in mHz - to get in Hz, divide by 1000
        [Pxx,~] = pmtm(ndat,nw,freq,1./ndt);
        H = dspdata.psd(Pxx,'Fs',1./ndt);
        psde = [H.Frequencies,H.Data];

        FA=zeros(nsamps,1);
        % Build FT of noise
        FA(1)=psde(1,2);
        FA(nsamps/2+1)=psde(end,2);
        FA(2:nsamps/2)=psde(2:end-1,2).*exp(1i*2*pi*random('unif',0,1,nsamps/2-1,1));
        FA(nsamps/2+2:end)=flipud(conj(FA(2:nsamps/2)));
        FA(:)=FA(:)/sqrt(mean(abs(FA(:)).^2)); % normalise to RMS=1
        FA(:)=FA(:)*sqrt(nsamps);
        %ifft to get noise data:
        rndat=ifft(FA);
        % get amplification factor to down-weight the noise  
        ampfac = maxab(avar.dataPSVSH(:,ic,ip))./rms(ndat)./maxab(idat); % N.B. rms of rndat=1
%         ampfac = 2*avar.SNR_stk(ic,ip)/maxab(trudata.(dtype)(ia).PSV(:,ic));
       
        % assemble final data
        odat(:,ic) = idat + noise_sta_deets.noiseup*rndat/ampfac;
        
        % posterior guess at SNR
%         inoise = trudata.(dtype)(ia).tt>=-10 & trudata.(dtype)(ia).tt<-3;
%         iSNR(ic) = max(odat(:,ic))./2./rms(odat(inoise,ic));

        
    end % loop on components
%          subplot(2,1,1),plot(trudata.(dtype)(ia).tt,trudata.(dtype)(ia).PSV)
%          subplot(2,1,2),plot(trudata.(dtype)(ia).tt,odat)
        
    trudata.(dtype)(ia).PSV = odat;
        
    end % loop on arrivals
end % loop on data types




end


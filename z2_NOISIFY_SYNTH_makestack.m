function [ trudata,par ] = z2_NOISIFY_SYNTH_makestack( trudata, par, noise_sta_deets )
%[ noisy_trudata,par ] = z2_NOISIFY_SYNTH( trudata, par, noise_sta_deets )
% 
% Function to add realistic noise to the synthetic data by grabbing an
% excerpt of noisy data from before the body wave arrivals (for Ps, Sp)
% from a real station, computing the power spectra, and adding noise to the
% synthetic data with this power spectrum (but randomised phase). 
% 
% This should be done before any of the data cleaning steps, as would be
% done for the real data. 

nwin = [-85 -25];
dwin = [-3 3];
proc_level = 4;

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


%% add realistic noise to the data
dtypes = fieldnames(trudata);
for id = 1:length(dtypes)
    dtype = dtypes{id};
    pdtyp = parse_dtype(dtype);
    if strcmp(pdtyp{1},'SW'), continue; end
    % make different data w/ different random noise - stops it being
    % dominated by single noisy peaks - simulates stack, effectively
    if strcmp(pdtyp{2},'Ps')
        trudata.(dtype)(2,1) = trudata.(dtype)(1);
    end
%     trudata.(dtype)(3,1) = trudata.(dtype)(1);
%     trudata.(dtype)(4,1) = trudata.(dtype)(1);
%     trudata.(dtype)(5,1) = trudata.(dtype)(1);
    
    for ia = 1:length(trudata.(dtype))

    % grab real data
    wd = pwd;
    cd(noise_sta_deets.datadir);
    datfiles = {};
    for gcwind = -4:4
        rdf = dir(sprintf('ngen_plotting_dat_%s_%s_%.0f_*.mat',noise_sta_deets.sta,noise_sta_deets.nwk,trudata.(dtype)(ia).gcarc+gcwind));
        datfiles = [datfiles;{rdf.name}']; %#ok<AGROW>
    end
    load(datfiles{1});
    cd(wd);

    ip = find(strcmp(pdtyp{2},plot_dat{1,3})); %#ok<*USENS>
    if strcmp(pdtyp{2},'Ps') 
        iparent = find(strcmp('P',plot_dat{2,3})); 
    elseif strcmp(pdtyp{2},'Sp')
        iparent = find(strcmp('SV',plot_dat{2,3})); 
    end


    odat = zeros(size(trudata.(dtype)(ia).PSV));
    for ic = 1:2
        idat = trudata.(dtype)(ia).PSV(:,ic);
        nsamps = length(idat);
            
            
        Ns = size(plot_dat{proc_level,ip}.data,2);
        psde_all = [];
        ampfac_all = [];
        for is = 1:Ns
            nind = plot_dat{proc_level,ip}.tt(:,is)>=nwin(1) & plot_dat{proc_level,ip}.tt(:,is)<nwin(2);
            dind_real = plot_dat{proc_level,ip}.tt(:,is)>=dwin(1) & plot_dat{proc_level,ip}.tt(:,is)<dwin(2);
            dind_syn = trudata.(dtype)(ia).tt>=dwin(1) & trudata.(dtype)(ia).tt<dwin(2);

                
            ndt = unique(round_level(diff(plot_dat{proc_level,ip}.tt(:,is)),0.0001));
            idt = 1./trudata.(dtype)(1).samprate;
            %% process noise series
            ndat_use = plot_dat{proc_level,ip}.data(nind,is,ic); 
            ndat_use = ndat_use(~isnan(ndat_use));
            ndat_use=detrend(ndat_use,'constant'); %detrend
%             ndat_use = plot_dat{proc_level,ip}.data(nind,is,1); % use pre-P noise
%             ndat_use = ndat_use(~isnan(ndat_use));
%             ndat_use=detrend(ndat_use,'constant'); %detrend
% 
%             ndat=ndat./sqrt(mean(ndat.^2)); %denorm

            %% calc individual noise spectra
            nw=4;
            freq=[1:floor(nsamps/2)+1]'./(nsamps*idt);%in mHz - to get in Hz, divide by 1000
            [Pxx,~] = pmtm(ndat_use,nw,freq,1./ndt);
%             Pxx = moving_average(Pxx,round(length(Pxx)/10)); % smooth power series!
            if strcmp(noise_sta_deets.noiseshape,'white')
                Pxx = mean(Pxx)*ones(size(Pxx));
            end

            H = dspdata.psd(Pxx,'Fs',1./ndt);
            psde_all(:,:,is) = [H.Frequencies,H.Data];
            % crucial line - the scaling of the noise
            ampfac_all(is) = max(abs(plot_dat{proc_level,ip}.data(dind_real,is,iparent)))./rms(ndat_use)./max(abs(trudata.(dtype)(ia).PSV(dind_syn,iparent)));
        end
        
        %% average the noise spectrum
        psde = [];
        psde(:,1) = psde_all(:,1,1);
        for iff = 10:length(psde)
            psde(iff,2) = mean_wtd(squeeze(psde_all(iff,2,:)),ampfac_all(:));
        end

        %% Make indiv random noise series
        ndat_all = zeros(length(idat),Ns);
        for is = 1:Ns        
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
%             rndat=detrend(rndat);
            % get amplification factor to down-weight the noise  
%             ampfac(ic) = mean(ampfac_all)/max(abs(idat)); % N.B. rms of rndat=1
    %         ampfac = 2*avar.SNR_stk(ic,ip)/maxab(trudata.(dtype)(ia).PSV(:,ic));
            ndat_all(:,is) = noise_sta_deets.noiseup*rndat/ampfac_all(is);
        end

        % assemble final data
        odat(:,ic) = idat + mean(ndat_all,2);
        
        % posterior guess at SNR
%         inoise = trudata.(dtype)(ia).tt>=-10 & trudata.(dtype)(ia).tt<-3;
%         iSNR(ic) = max(odat(:,ic))./2./rms(odat(inoise,ic));

        
    end % loop on components

%     figure(4)
%     subplot(2,1,1),plot(trudata.(dtype)(ia).tt,trudata.(dtype)(ia).PSV), ylim(max(max(abs(odat)))*[-1 1])
%     subplot(2,1,2),plot(trudata.(dtype)(ia).tt,odat), ylim(max(max(abs(odat)))*[-1 1])
%         
    trudata.(dtype)(ia).PSV = odat;
        
    end % loop on arrivals
end % loop on data types




end


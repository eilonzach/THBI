function [ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,data,datahparm,par )
%  [ log_likelihood,misfit ] = b5_CALC_LIKELIHOOD( misfit,data,datahparm,par )
% 
%   Function to calculate the likelihood function (the inverse exponential
%   of the error-weighted misfit function)
% 
% INPUTS:
%  misfit = misfit structure with squared misfit of each data type
%  data   = true data structure, including estimated sigma of each data
%           type
%  par    = input parameters, including what data types to solve for

% FOR each datatype, do the following:
%  1) get data error
%  2) Number of datapoints
%  3) chi2 of each data type
%  4) save chi2 into misfit structure
%  5) save rms into misfit structure
%  6) calc log_likelihood
% 	   ==> likelihood of each data type
%         sq2p = sqrt(2*pi);
%             we will ignore a constant factor of 
%             -0.5*N_ps*log(sq2p) -0.5*N_sp*log(sq2p) -0.5*N_SW*log(sq2p) 
%             because it never changes, and so does not affect results

%% SW
if any(strcmp(par.inv.datatypes,'SW'))
    sigma_SW = datahparm.sigmaSW;
    N_SW = length(data.SW.periods);
    chi2_SW = misfit.SW/sigma_SW/sigma_SW; 
    misfit.chi2_SW = chi2_SW;
    misfit.rms_SW = sqrt(misfit.SW/N_SW);
    log_likelihood_SW = -N_SW*log(sigma_SW) - 0.5*chi2_SW;
else
    chi2_SW = 0;
    log_likelihood_SW = 0; 
end

%% Ps RF
if any(strcmp(par.inv.datatypes,'PsRF'))
    sigma_ps = datahparm.sigmaPsRF;
    % N_ps = data.PsRF.nsamp;
    M_ps = length(data.PsRF);
    N_ps = zeros(M_ps,1);
    for itr = 1:length(data.PsRF)
        N_ps(itr) = mean([scdofcalc(data.PsRF(itr).ZRT(:,1)),scdofcalc(data.PsRF(itr).ZRT(:,2))]);
    end
    chi2_ps = misfit.PsRF/sigma_ps/sigma_ps;
    misfit.N_ps = N_ps;
    misfit.chi2_ps = chi2_ps;
    misfit.rms_ps = sqrt(misfit.PsRF./N_ps);
    log_likelihood_ps = sum(-N_ps*log(sigma_ps) - 0.5*chi2_ps)./M_ps;
else
    M_ps = 1;
    chi2_ps = 0;
    log_likelihood_ps = 0; 
end

%% Sp RF
if any(strcmp(par.inv.datatypes,'SpRF')), 
    sigma_sp = datahparm.sigmaSpRF;
    % N_sp = data.SpRF.nsamp;
    M_sp = length(data.SpRF);
    N_sp = zeros(M_sp,1);
    for itr = 1:length(data.SpRF)
        N_sp(itr) = mean([scdofcalc(data.SpRF(itr).ZRT(:,1)),scdofcalc(data.SpRF(itr).ZRT(:,2))]);
    end
    chi2_sp = misfit.SpRF/sigma_sp/sigma_sp;
	misfit.N_sp = N_sp;
    misfit.chi2_sp = chi2_sp;
    misfit.rms_sp = sqrt(misfit.SpRF./N_sp);
    log_likelihood_sp = sum(-N_sp*log(sigma_sp) - 0.5*chi2_sp)./M_sp;
else
    M_sp = 1;
    chi2_sp = 0;
    log_likelihood_sp = 0; 
end

%% Ps RF - low_f
if any(strcmp(par.inv.datatypes,'PsRF_lo'))
    sigma_ps_lo = datahparm.sigmaPsRF_lo;
    % N_ps_lo = data.PsRF_lo.nsamp;
    M_ps_lo = length(data.PsRF_lo);
    N_ps_lo = zeros(M_ps_lo,1);
    for itr = 1:length(data.PsRF_lo)
        N_ps_lo(itr) = mean([scdofcalc(data.PsRF_lo(itr).ZRT(:,1)),scdofcalc(data.PsRF_lo(itr).ZRT(:,2))]);
    end
    chi2_ps_lo = misfit.PsRF_lo/sigma_ps_lo/sigma_ps_lo;
    misfit.N_ps_lo = N_ps_lo;
    misfit.chi2_ps_lo = chi2_ps_lo;
    misfit.rms_ps_lo = sqrt(misfit.PsRF_lo./N_ps_lo);
    log_likelihood_ps_lo = sum(-N_ps_lo*log(sigma_ps_lo) - 0.5*chi2_ps_lo)./M_ps_lo;
else
    M_ps_lo = 1;
    chi2_ps_lo = 0;
    log_likelihood_ps_lo = 0;
end

%% Sp RF - low_f
if any(strcmp(par.inv.datatypes,'SpRF_lo'))
    sigma_sp_lo = datahparm.sigmaSpRF_lo;
    % N_sp_lo = data.SpRF_lo.nsamp;
    M_sp_lo = length(data.SpRF_lo);
    N_sp_lo = zeros(M_sp_lo,1);
    for itr = 1:length(data.SpRF_lo)
        N_sp_lo(itr) = mean([scdofcalc(data.SpRF_lo(itr).ZRT(:,1)),scdofcalc(data.SpRF_lo(itr).ZRT(:,2))]);
    end
    chi2_sp_lo = misfit.SpRF_lo/sigma_sp_lo/sigma_sp_lo;
    misfit.N_sp_lo = N_sp_lo;
    misfit.chi2_sp_lo = chi2_sp_lo;
    misfit.rms_sp_lo = sqrt(misfit.SpRF_lo./N_sp_lo);
    log_likelihood_sp_lo = sum(-N_sp_lo*log(sigma_sp_lo) - 0.5*chi2_sp_lo)./M_sp_lo;
else
    M_sp_lo = 1;
    chi2_sp_lo = 0;
    log_likelihood_sp_lo = 0; 
end


%% save overall misfit
misfit.chi2 = sum(chi2_SW) + sum(chi2_ps)./M_ps + sum(chi2_sp)./M_sp + sum(chi2_ps_lo)./M_ps_lo + sum(chi2_sp_lo)./M_sp_lo;
misfit = orderfields(misfit,sort(fieldnames(misfit))); % order the fields for neatness

%% calc. overall log likelihood
log_likelihood = log_likelihood_ps + log_likelihood_sp + log_likelihood_SW + log_likelihood_ps_lo + log_likelihood_sp_lo;


end


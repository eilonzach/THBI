function [ noisy_trudata,par ] = z2_NOISIFY_SYNTH( trudata, par, noise_sta_deets )
%[ noisy_trudata,par ] = z2_NOISIFY_SYNTH( trudata, par, noise_sta_deets )
% 
% Function to add realistic noise to the synthetic data by grabbing an
% excerpt of noisy data from before the body wave arrivals (for Ps, Sp)
% from a real station, computing the power spectra, and adding noise to the
% synthetic data with this power spectrum (but randomised phase). 
% 
% This should be done before any of the data cleaning steps, as would be
% done for the real data. 

if nargin < 3 || isempty(noise_sta_deets)
noise_sta_deets = struct('projname','WYOMING','sta','RSSD','nwk','IU','gc','70');
end

noisedata = load_data(noise_sta_deets.projname,noise_sta_deets.sta,noise_sta_deets.nwk,noise_sta_deets.gc);






end


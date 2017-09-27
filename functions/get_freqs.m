function [ freqs, periods ] = get_freqs(datadir)
%[ freqs, periods ] = get_freqs(datadir)
%   Quick function to spit out all of the frequencies available. N.B.
%   assumes files are in the format "gridvalue.velFFFkern.dat", where "FFF"
%   is a three-character value for the frequency (e.g. "006") in mHz.

if nargin<1 || isempty(datadir)
    datadir = '~/Dropbox/Dave_Li_phV/2D_phase_velocities/'; % need final slash
end

phVfiles = dir(datadir); 
freqs = nan(length(phVfiles),1);
for ii = 1:length(phVfiles), 
    try freqs(ii) = eval(phVfiles(ii).name(14:16)); end, 
end;
freqs(isnan(freqs))=[]; % frequencies are in mHz

freqs = freqs*0.001;
periods = 1./freqs;


end


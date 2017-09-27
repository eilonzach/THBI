function [ phV, periods ] = RD_1DphV
%[ phV, periods ] = RD_1DphV
%   Function to grab the 1D phase velocity profile from a digitised version
%   of the blue curve in Riddhi Dave and Aibing Li's 2016 paper's online
%   supplemental figure DR3_A


a = load('~/Dropbox/Dave_Li_phV/Wyoming_phV.txt','-ascii');

freqs = round_level(1./a(:,1),0.001);
periods = 1./freqs;

phV = a(:,2);



end


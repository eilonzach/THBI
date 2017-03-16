function [ dat_out ] = predat_process( dat_in,cp,ifnorm,ifdec )
%[ dat_out ] = predat_process( dat_in,cp,ifnorm,ifdec )
%   
%  Function to process predicted data by cutting, tapering, and normalising
%  to unit energy. Then decimates to 2*Nyquist.
% 
% INPUT:
%  dat_in  - structure with fields:
%             .ZRT - data in columns, ordered Z,R,T
%             .tt  - tt vector, w/ significant padding before phase arrival
%             .nsamp - number of samples
%  cp      - cleaning_parm structure, as for use with data_clean
%  ifnorm  - option [true/false] to normalise data to unit energy
%  ifdec   - option [true/false] to decimate data to 2*fNyq
%       
% OUTPUT:
%  dat_out - similar structure, but with the data having been cut,
%               tapered with a 5 second hanning window at either end, and
%               normalised to unit energy.
% 
%  Z. Eilon 09/2016

if nargin<3 || isempty(ifnorm)
    ifnorm=false;
end
if nargin<4 || isempty(ifdec)
    ifdec=false;
end

dat_out = dat_in;

cp.samprate = dat_in.samprate;

%% clean, filter, taper
dat_out.ZRT = data_clean(dat_out.ZRT,cp);

%% normalise to unit energy
if ifnorm && ~isempty(dat_out.ZRT)
    normf_ps = dat_out.ZRT(:,1)'*dat_out.ZRT(:,1) + ...
               dat_out.ZRT(:,2)'*dat_out.ZRT(:,2) + ...
               dat_out.ZRT(:,3)'*dat_out.ZRT(:,3);
    dat_out.ZRT = dat_out.ZRT/sqrt(normf_ps);
end
    
%% decimate
if ifdec && ~isempty(dat_out.ZRT)
    resamprate = cp.fhi*4;
    tt_new = [dat_out.tt(1):1./resamprate:dat_out.tt(end)]';
    dat_out.ZRT = interp1(dat_out.tt,dat_out.ZRT,tt_new);
    dat_out.tt = tt_new;
    dat_out.nsamp = length(dat_out.tt);
end

end

